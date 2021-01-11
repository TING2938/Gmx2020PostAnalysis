
/* run as main() for gromacs-2018.* code */
#define gmx_main(func)										\
		int func(int argc, char** argv);					\
		int main(int argc, char** argv)						\
		{													\
			return gmx_run_cmain(argc, argv, &func);		\
		}													\
		int func(int argc, char** argv)

/* function get current dir */
#include <stdio.h>  /* defines FILENAME_MAX */
#ifdef _WIN32
	#include <direct.h>
	#define GetCurrentDir _getcwd
#else
	#include <unistd.h>
	#define GetCurrentDir getcwd
#endif