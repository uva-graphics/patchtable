
/* Read/write 1 and 3 channel PFM files, public domain Connelly Barnes 2007. */

#ifndef _pfm_h
#define _pfm_h

int is_little_endian();
float *read_pfm_file(const char *filename, int *w, int *h);
void write_pfm_file(const char *filename, float *depth, int w, int h);
float *read_pfm_file3(const char *filename, int *w, int *h);
void write_pfm_file3(const char *filename, float *depth, int w, int h);

#endif
