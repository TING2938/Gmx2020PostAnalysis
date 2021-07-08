#include <itp/gmx>

#define LOG(a) printf("%s: %d", #a, a)

gmx_main(func)
{
    int fd;
    LOG(fd);

    return 0;
}