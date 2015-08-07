
#ifndef _product_quantize_h
#define _product_quantize_h

#include <opencv2/core/core.hpp>

#define TABLE_PRODUCT_QUANTIZE_VERBOSE    0             /* Extra verbosity for debugging (should be turned off for checked in code) */

template<class real, class itype=TABLE_DEFAULT_ITYPE>
class ProductQuantizer { public:
    PatchTableParams *p;
    const Array<real> &wh0;
    const vector<int> &nslices;
    const Array<itype> *allowed_patches;
    
    int product_count;
    vector<int> dim_lo, dim_hi, dim_count;
    vector<int> cluster_count;
    
    vector<shared_ptr<Array<real> > > cluster_centers;     /* Cluster center matrix, clusters x ndims, length is product_count. */
    vector<shared_ptr<Array<real> > > pairwise_distance;   /* Cluster pairwise distance matrix, length is product_count. */
    vector<shared_ptr<Array<int> > > cluster_knn;          /* clusters x k nearest clusters for each cluster, product_count length */
    vector<shared_ptr<Array<int> > > grid_to_cluster;      /* Map from grid cell to cluster index, length is product_count */
    vector<vector<real> > grid_min;                        /* Min coord for any cluster */
    vector<vector<real> > grid_max;                        /* Max coord for any cluster */
    
    vector<int> temp_quantize_index;
    
    ProductQuantizer(PatchTableParams *p_, const Array<real> &wh0_, const vector<int> &nslices_, const Array<itype> *allowed_patches_=NULL) :p(p_), wh0(wh0_), nslices(nslices_), allowed_patches(allowed_patches_) {
        double T0_constructor = wall_time();
        product_count = DIV_ROUND_UP(wh0.channels(), p->product_quantize_dims);
        
        ASSERT2(wh0.channels() == nslices.size(), "ProductQuantizer: expected wh0.channels() == nslices.size()");
        
        temp_quantize_index.resize(product_count);
        
        for (int i = 0; i < product_count; i++) {
            int lo = i*p->product_quantize_dims;
            dim_lo.push_back(lo);
            int hi = (i+1)*p->product_quantize_dims;
            if (hi > wh0.channels()) { hi = wh0.channels(); }
            dim_hi.push_back(hi);
            dim_count.push_back(hi-lo);
            
            int current_cluster_count = 1;
            for (int j = dim_lo[i]; j < dim_hi[i]; j++) {
                current_cluster_count *= nslices[j];
            }
            cluster_count.push_back(current_cluster_count);
        }
        
        int point_count = 0;
        TABLE_FOR_ALLOWED_PATCHES() {
            point_count++;
        } TABLE_END_FOR_ALLOWED_PATCHES();
        
        double T_kmeans = 0;
        double T_pairwise_distance = 0;
        double T_grid_to_cluster = 0;
        
        grid_min.resize(product_count);
        grid_max.resize(product_count);
        
        for (int i = 0; i < product_count; i++) {
            int lo = dim_lo[i];
            int hi = dim_hi[i];
            int dims = dim_count[i];
            cv::Mat points(point_count, dims, CV_32FC1);
            
            int current_point = 0;
            TABLE_FOR_ALLOWED_PATCHES() {
                float *row = points.ptr<float>(current_point);
                real *wh0_p = &wh0.get_nearest(y, x, lo);
                for (int j = 0; j < dims; j++) {
                    row[j] = wh0_p[j];
                }
                current_point++;
            } TABLE_END_FOR_ALLOWED_PATCHES();
            ASSERT2(current_point == point_count, "expected current_point == point_count");
            
            int cluster_count_i = cluster_count[i];
            
            cv::Mat labels, centers;
            cv::TermCriteria criteria(cv::TermCriteria::EPS+cv::TermCriteria::COUNT, p->kmeans_max_iters, p->kmeans_eps);
            double T0_kmeans = wall_time();
            kmeans(points, cluster_count_i, labels, criteria, p->kmeans_attempts, cv::KMEANS_PP_CENTERS, centers);
            T_kmeans += wall_time() - T0_kmeans;
            ASSERT2(labels.rows == point_count, "expected labels.rows == point_count");
            ASSERT2(labels.cols == 1, "expected labels.cols == 1");
            ASSERT2(centers.rows == cluster_count_i, "expected clusters.rows == cluster_count_i");
            ASSERT2(centers.cols == dims, "expected centers.cols == dims");
            
            grid_min[i].resize(dims);
            grid_max[i].resize(dims);
            vector<real> &grid_min_v(grid_min[i]);
            vector<real> &grid_max_v(grid_max[i]);
            for (int j = 0; j < dims; j++) {
                grid_min_v[j] = 1e100;
                grid_max_v[j] = -1e100;
            }
            
            shared_ptr<Array<real> > centers_ours(make_shared<Array<real> >(centers.rows, centers.cols));
            for (int row = 0; row < centers.rows; row++) {
                float *centers_row = centers.ptr<float>(row);
                real *centers_ours_row = &centers_ours->get_nearest(row, 0);
                for (int col = 0; col < centers.cols; col++) {
                    float value = centers_row[col];
                    if (value < grid_min_v[col]) { grid_min_v[col] = value; }
                    if (value > grid_max_v[col]) { grid_max_v[col] = value; }
                    centers_ours_row[col] = value;
                }
            }
            cluster_centers.push_back(centers_ours);
            
            double T0_pairwise_distance = wall_time();
            shared_ptr<Array<real> > D(make_shared<Array<real> >(cluster_count_i, cluster_count_i));
            for (int j = 0; j < cluster_count_i; j++) {
                for (int k = 0; k < cluster_count_i; k++) {
                    real *centers_ours_j = &centers_ours->get_nearest(j, 0);
                    real *centers_ours_k = &centers_ours->get_nearest(k, 0);
                    real dist = 0;
                    for (int l = 0; l < dims; l++) {
                        real delta = centers_ours_j[l] - centers_ours_k[l];
                        dist += delta*delta;
                    }
                    (*D)(j, k) = dist;
                }
            }
            pairwise_distance.push_back(D);
            
            int knn = MIN(p->product_quantize_knn, cluster_count_i-1);
            shared_ptr<Array<int> > knn_arr(make_shared<Array<int> >(cluster_count_i, knn));
            vector<pair<real, int> > D_row(cluster_count_i-1);
            for (int j = 0; j < cluster_count_i; j++) {
                for (int k = 0; k < j; k++) {
                    D_row[k] = pair<real, int>((*D)(j, k), k);
                }
                for (int k = j+1; k < cluster_count_i; k++) {
                    D_row[k-1] = pair<real, int>((*D)(j, k), k);
                }

                std::partial_sort(D_row.begin(), D_row.begin() + knn, D_row.end());
                for (int k = 0; k < knn; k++) {
                    (*knn_arr)(j, k) = D_row[k].second;
                }
            }
            cluster_knn.push_back(knn_arr);
            
            T_pairwise_distance += wall_time() - T0_pairwise_distance;
            
            if (p->product_quantize_log) {
                char buf[256];
                sprintf(buf, "product_quantize_%d_centers.txt", i);
                FILE *f = fopen(buf, "wt");
                
                for (int row = 0; row < centers.rows; row++) {
                    for (int col = 0; col < centers.cols; col++) {
                        fprintf(f, "%f ", double(centers_ours->get_nearest(row, col)));
                    }
                    fprintf(f, "\n");
                }
                
                fclose(f);
                
                sprintf(buf, "product_quantize_%d_points.txt", i);
                f = fopen(buf, "wt");
                
                for (int j = 0; j < point_count; j++) {
                    float *row = points.ptr<float>(j);
                    for (int k = 0; k < dims; k++) {
                        fprintf(f, "%f ", double(row[k]));
                    }
                    fprintf(f, "\n");
                }
                
                fclose(f);
                
                sprintf(buf, "product_quantize_%d_labels.txt", i);
                f = fopen(buf, "wt");
                
                for (int j = 0; j < point_count; j++) {
                    fprintf(f, "%d\n", labels.at<int>(j));
                }
                
                fclose(f);

                sprintf(buf, "product_quantize_%d_dist.txt", i);
                f = fopen(buf, "wt");
                
                for (int j = 0; j < cluster_count_i; j++) {
                    for (int k = 0; k < cluster_count_i; k++) {
                        fprintf(f, "%f ", (*D)(j, k));
                    }
                    fprintf(f, "\n");
                }
                
                fclose(f);

                sprintf(buf, "product_quantize_%d_knn.txt", i);
                f = fopen(buf, "wt");
                
                for (int j = 0; j < cluster_count_i; j++) {
                    for (int k = 0; k < knn; k++) {
                        fprintf(f, "%d ", (*knn_arr)(j, k));
                    }
                    fprintf(f, "\n");
                }
                
                fclose(f);

            }
            
            double T0_grid_to_cluster = wall_time();
            shared_ptr<Array<int> > G_ptr(make_shared<Array<int> >());
            vector<int> G_sizes(dims, p->product_quantize_mapn);
            Array<int> &G(*G_ptr.get());
            G.resize(G_sizes);
            
            flann::Matrix<real> flann_cluster_centers(centers_ours->data, centers_ours->height(), centers_ours->width());
            flann::Index<flann::L2<real> > flann_index(flann_cluster_centers, flann::KDTreeSingleIndexParams());
            flann_index.buildIndex();
            
            vector<real> query_vector(dims);
            
            vector<vector<int> > matched_indices_mat;
            vector<vector<real> > matched_dists_mat;
            flann::Matrix<real> query_matrix(&query_vector[0], 1, query_vector.size());
            matched_indices_mat.resize(1);
            matched_indices_mat[0].resize(1);
            matched_dists_mat.resize(1);
            matched_dists_mat[0].resize(1);
            flann::SearchParams search_params(p->flann_checks);
            search_params.eps = 0.0;
            
            for (int j = 0; j < G.nelems; j++) {
                for (int l = 0; l < dims; l++) {
                    int j_dim = (j / G.stride[l]) % G.sizes[l];
                    real t = (j_dim+0.5) / G.sizes[l];
                    query_vector[l] = grid_min_v[l] + (grid_max_v[l]-grid_min_v[l]) * t;
                }
                flann_index.knnSearch(query_matrix, matched_indices_mat, matched_dists_mat, 1, search_params);
                int idx = matched_indices_mat[0][0];
                ASSERT2(in_bounds(idx, cluster_count_i), "expected grid_to_cluster matched index to be in cluster_count_i bounds");
                G.data[j] = idx;
            }

            if (p->product_quantize_log) {
                char buf[256];
                sprintf(buf, "product_quantize_%d_grid_to_cluster.txt", i);
                FILE *f = fopen(buf, "wt");
                
                for (int j = 0; j < G.nelems; j++) {
                    fprintf(f, "%d ", G.data[j]);
                }
                
                fclose(f);
            }
            
            grid_to_cluster.push_back(G_ptr);
            T_grid_to_cluster += wall_time()-T0_grid_to_cluster;
            
            if (p->verbose) {
                printf("product i=%d, T_kmeans: %f, T_pairwise_distance: %f, T_grid_to_cluster: %f\n", i, T_kmeans, T_pairwise_distance, T_grid_to_cluster);
            }
        }
        if (p->verbose) {
            printf("ProductQuantizer: %f secs\n", wall_time()-T0_constructor);
        }
    }

    /* Map real vector into int indices of length product_count */
    void quantize(real *v, int *index) {
        int ans = 0;
        for (int i = 0; i < product_count; i++) {
            int j = 0;
            int dims = dim_count[i];
            int lo = dim_lo[i];
            vector<real> &grid_min_v(grid_min[i]);
            vector<real> &grid_max_v(grid_max[i]);
            Array<int> &G = *grid_to_cluster[i].get();
            
            for (int l = 0; l < dims; l++) {
                real t = (v[lo+l] - grid_min_v[l]) / (grid_max_v[l] - grid_min_v[l]);
                if (t < 0) { t = 0; }
                else if (t >= 1) { t = 1 - 1e-7; }
                int j_dim = int(G.sizes[l]*t);
                ASSERT2(in_bounds(j_dim, G.sizes[l]), "expected j_dim in bounds in quantize");
                j += j_dim * G.stride[l];
            }
            ASSERT2(in_bounds(j, G.nelems), "expected j in bounds in quantize");
            
            index[i] = G.data[j];
        }
    }

    /* Map real vector into single int index */
    int quantize(real *v) {
#if TABLE_PRODUCT_QUANTIZE_VERBOSE
        if (p->verbose >= 2) {
            printf("beginning of quantize, product_count=%d\n", product_count);
        }
#endif
        int ans = 0;
        int stride = 1;
        
        quantize(v, &temp_quantize_index[0]);
#if TABLE_PRODUCT_QUANTIZE_VERBOSE
        if (p->verbose >= 2) {
            printf("called quantize(real *, int *), calculating a single int index\n");
        }
#endif
        for (int i = product_count-1; i >= 0; i--) {
            ans += temp_quantize_index[i] * stride;
            stride *= cluster_count[i];
        }
#if TABLE_PRODUCT_QUANTIZE_VERBOSE
        if (p->verbose >= 2) {
            printf("nslices: %s\n", vector_to_str_int(nslices).c_str());
            printf("cluster_count: %s\n", vector_to_str_int(cluster_count).c_str());
            
            printf("temp_quantize_index: ");
            for (int i = 0; i < product_count; i++) {
                printf("%d ", temp_quantize_index[i]);
            }
            printf("product of all dims: %d\n", stride);
            
            stride = 1;
            printf("stride products in reverse order: ");
            for (int i = product_count-1; i >= 0; i--) {
                printf("%d\n", stride);
                stride *= cluster_count[i];
            }
        }
#endif

        return ans;
    }
    
    /* Map single int index back to cluster center in product space (with length wh0.channels()) */
    void get_cluster_center(int index, real *center) {
        int stride = 1;
        
        for (int i = product_count-1; i >= 0; i--) {
            int cluster_index = (index/stride) % cluster_count[i];
            stride *= cluster_count[i];

            int lo = dim_lo[i];
            int dims = dim_count[i];
            Array<real> &C(*cluster_centers[i]);
            ASSERT2(in_bounds(cluster_index, C.height()), "expected cluster_index in bounds");
            real *C_row = &C(cluster_index, 0);
            for (int j = 0; j < dims; j++) {
                center[j+lo] = C_row[j];
            }
        }
    }
};

#endif
