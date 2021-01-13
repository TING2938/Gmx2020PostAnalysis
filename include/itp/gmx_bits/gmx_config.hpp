
/* run as main() for gromacs-2018.* code */
#define gmx_main(func)										\
		int func(int argc, char** argv);					\
		int main(int argc, char** argv)						\
		{													\
			return gmx_run_cmain(argc, argv, &func);		\
		}													\
		int func(int argc, char** argv)

