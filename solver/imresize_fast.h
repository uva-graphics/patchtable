template<int n>
double lanczos(double x) {
    if (x >= -n && x <= n) {
        return sinc(x) * sinc(x/double(n));
    }
    return 0.0;
}

#define IMRESIZE_NLANCZOS 1024

template<class real, int channels=3, int kernel_n=4>
void imresize_fast(const Array<real> &input, Array<real> &output, int w, int h, bool prefilter_gaussian=true) {
    ASSERT2(input.channels() == channels, "expected input.channels() matching template channel parameter");
    int in_w = input.width(), in_h = input.height();
    double T0 = wall_time();
    
    Array<real> temp;
    temp.assign(input);
    bool is_prefilter = false;
    if (prefilter_gaussian && (input.height() > h || input.width() > w)) {
        is_prefilter = true;
        double sigma_y = double(input.height())/double(h);
        double sigma_x = double(input.width())/double(w);
        if (sigma_y < 1) { sigma_y = 0; }
        if (sigma_x < 1) { sigma_x = 0; }
        sigma_x *= 0.5;
        sigma_y *= 0.5;
        double sigma_both = (sigma_x+sigma_y)*0.5;
        gaussian_blur<real, -1, channels>(temp, sigma_both);
//        save_color_image<real>(temp, "temp_blurred_iir.png");
        //cv::GaussianBlur(input, input, cv::Size(0, 0), sigma_x, sigma_y, cv::BORDER_REPLICATE);
    }
    printf("imresize_fast: prefilter: %f secs\n", wall_time()-T0);
    
    const Array<real> &temp_array = is_prefilter ? temp: input;
    output.resize(h, w, channels);

    static double lanczos_discretized[IMRESIZE_NLANCZOS];
    static bool lanczos_init = false;
    if (!lanczos_init) {
        lanczos_init = true;
        for (int i = 0; i < IMRESIZE_NLANCZOS; i++) {
            double v_arg = kernel_n * i / (IMRESIZE_NLANCZOS-1);
            lanczos_discretized[i] = (i == (IMRESIZE_NLANCZOS - 1) ? 0.0: lanczos<kernel_n>(v_arg));
        }
    }

    #pragma omp parallel for
    for (int y = 0; y < h; y++) {
        double yfloat = y*1.0*(in_h-1)/(h-1);
        int yi(yfloat);
        double yf(yfloat-yi);
        
        double kernel_weight_y[kernel_n*2+2];
        double kernel_weight_x[kernel_n*2+2];

        double kernel_weight_ysum = 0.0;
        for (int dy = -kernel_n-1; dy <= kernel_n; dy++) {
            double yarg = dy - yf;
            if (yarg < 0) { yarg = -yarg; }
            if (yarg > kernel_n) { yarg = kernel_n; }
            int y_index = yarg * (IMRESIZE_NLANCZOS-1) / kernel_n;
            double kernel_weight_y_current = lanczos_discretized[y_index];
            kernel_weight_ysum += kernel_weight_y_current;
            kernel_weight_y[dy+kernel_n+1] = kernel_weight_y_current;
        }
        
        for (int x = 0; x < w; x++) {
            double xfloat = x*1.0*(in_w-1)/(w-1);
            int xi(xfloat);
            double xf(xfloat-xi);
            
            double kernel_weight_xsum = 0.0;
            
            for (int dx = -kernel_n-1; dx <= kernel_n; dx++) {
                double xarg = dx - xf;
                if (xarg < 0) { xarg = -xarg; }
                if (xarg > kernel_n) { xarg = kernel_n; }
                int x_index = xarg * (IMRESIZE_NLANCZOS-1) / kernel_n;
                double kernel_weight_x_current = lanczos_discretized[x_index];
                kernel_weight_xsum += kernel_weight_x_current;
                kernel_weight_x[dx+kernel_n+1] = kernel_weight_x_current;
            }
            
            double c_out[channels] = { 0.0 };
            for (int dy = -kernel_n-1; dy <= kernel_n; dy++) {
                double kernel_weight_y_current = kernel_weight_y[dy+kernel_n+1];
                
                int yp = yi + dy;
                if (yp < 0) { yp = 0; }
                else if (yp >= in_h) { yp = in_h-1; }
                real *temp_row = temp_array.data + yp * temp_array.stride[0];
                
                for (int dx = -kernel_n-1; dx <= kernel_n; dx++) {
                    double kernel_weight = kernel_weight_y_current * kernel_weight_x[dx+kernel_n+1];
                    //double kernel_weight = kernel_weight_y * lanczos<4>(xarg);
                    
                    int xp = xi + dx;
                    if (xp < 0) { xp = 0; }
                    else if (xp >= in_w) { xp = in_w-1; }
                    
                    real *temp_pixel = temp_row + xp * temp_array.stride[1];
                    
                    for (int c = 0; c < channels; c++) {
                        c_out[c] += kernel_weight * temp_pixel[c];
                    }
                }
            }

            double c_scale = 1.0/(kernel_weight_xsum*kernel_weight_ysum);
            for (int c = 0; c < channels; c++) {
                output.get_nearest(y, x, c) = c_out[c] * c_scale;
            }
        
        }
    }
    //cv::resize(input, output, cv::Size(w, h), 0.0, 0.0, interpolation);
    //return Array<real>(output);
}

