#include <sys/time.h>
#include <sys/resource.h>
#include "header.h"

long double get_time()
{
	struct rusage buf;
	
	getrusage (RUSAGE_SELF, &buf);
	return (double)buf.ru_utime.tv_sec + (double)buf.ru_utime.tv_usec / 1000000.0;
}

long double get_full_time()
{
    struct timeval buf;
    gettimeofday(&buf, 0);
    return (double)buf.tv_sec + (double)buf.tv_usec/1000000.0;
}
