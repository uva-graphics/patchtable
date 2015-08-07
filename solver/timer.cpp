//#include "array.h"

#if VISUAL_STUDIO_WORKAROUND

#include <windows.h>
bool wall_time_initialized = false;
LARGE_INTEGER wall_time_rate, wall_time_start;

double wall_time() {
    if (!wall_time_initialized) {
        wall_time_initialized = true;
		QueryPerformanceFrequency(&wall_time_rate);
        //ASSERT2(QueryPerformanceFrequency(&wall_time_rate), "timer did not initialize");
        QueryPerformanceCounter(&wall_time_start);
    }
    LARGE_INTEGER t;
    QueryPerformanceCounter(&t);
    return (t.QuadPart - wall_time_start.QuadPart) * 1.0 / wall_time_rate.QuadPart;
}

#else
#include <sys/time.h>
#include <cstddef>

double wall_time() {
    timeval t;
    gettimeofday(&t, NULL);
    return t.tv_sec + t.tv_usec*1e-6;
}
#endif
