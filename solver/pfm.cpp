
/* Read/write 1 and 3 channel PFM files, public domain Connelly Barnes 2007. */

#include <stdlib.h>
#include <string.h>
#include <stdio.h>

typedef unsigned char byte;

int is_little_endian() {
  if (sizeof(float) != 4) { printf("Bad float size.\n"); exit(1); }
  byte b[4] = { 255, 0, 0, 0 };
  return *((float *) b) < 1.0;
}

float *read_pfm_file3(const char *filename, int *w, int *h) {
  char buf[256];
  FILE *f = fopen(filename, "rb");
  fscanf(f, "%s\n", buf);
  if (strcmp(buf, "PF") != 0) {
    printf("Not a 3 channel PFM file.\n"); return NULL;
  }
  fscanf(f, "%d %d\n", w, h);
  double scale = 1.0;
  fscanf(f, "%lf\n", &scale);
  int little_endian = 0;
  if (scale < 0.0) {
    little_endian = 1;
    scale = -scale;
  }
  int channels = 3;
  byte *data = new byte[(*w)*(*h)*4*channels];
  float *depth = new float[(*w)*(*h)*channels];
  int count = fread((void *) data, 4, (*w)*(*h)*channels, f);
  if (count != (*w)*(*h)*channels) {
    printf("Error reading PFM file %d %d %d %d %d.\n", count, *w, *h, channels, (*w)*(*h)*channels); return NULL;
  }
  int native_little_endian = is_little_endian();
  for (int i = 0; i < (*w)*(*h)*channels; i++) {
    byte *p = &data[i*4];
    if (little_endian != native_little_endian) { 
      byte temp;
      temp = p[0]; p[0] = p[3]; p[3] = temp;
      temp = p[1]; p[1] = p[2]; p[2] = temp;
    }
    depth[i] = *((float *) p);
  }
  fclose(f);
  delete[] data;
  return depth;
}

void write_pfm_file3(const char *filename, float *depth, int w, int h) {
  FILE *f = fopen(filename, "wb");

  double scale = 1.0;
  if (is_little_endian()) { scale = -scale; }

  fprintf(f, "PF\n%d %d\n%lf\n", w, h, scale);
  int channels = 3;
  for (int i = 0; i < w*h*channels; i++) {
    float d = depth[i];
    fwrite((void *) &d, 1, 4, f);
  }
  fclose(f);
}

float *read_pfm_file(const char *filename, int *w, int *h) {
  char buf[256];
  FILE *f = fopen(filename, "rb");
  fscanf(f, "%s\n", buf);
  if (strcmp(buf, "Pf") != 0) {
    //printf("Not a 1 channel PFM file.\n");
    return NULL;
  }
  fscanf(f, "%d %d\n", w, h);
  double scale = 1.0;
  fscanf(f, "%lf\n", &scale);
  int little_endian = 0;
  if (scale < 0.0) {
    little_endian = 1;
    scale = -scale;
  }
  byte *data = new byte[(*w)*(*h)*4];
  float *depth = new float[(*w)*(*h)];
  int count = fread((void *) data, 4, (*w)*(*h), f);
  if (count != (*w)*(*h)) {
    printf("Error reading PFM file.\n"); return NULL;
  }
  int native_little_endian = is_little_endian();
  for (int i = 0; i < (*w)*(*h); i++) {
    byte *p = &data[i*4];
    if (little_endian != native_little_endian) { 
      byte temp;
      temp = p[0]; p[0] = p[3]; p[3] = temp;
      temp = p[1]; p[1] = p[2]; p[2] = temp;
    }
    depth[i] = *((float *) p);
  }
  fclose(f);
  delete[] data;
  return depth;
}

void write_pfm_file(const char *filename, float *depth, int w, int h) {
  FILE *f = fopen(filename, "wb");

  double scale = 1.0;
  if (is_little_endian()) { scale = -scale; }

  fprintf(f, "Pf\n%d %d\n%lf\n", w, h, scale);
  for (int i = 0; i < w*h; i++) {
    float d = depth[i];
    fwrite((void *) &d, 1, 4, f);
  }
  fclose(f);
}
