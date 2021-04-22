
namespace itp
{
	/**
	 * @brief Gromacs轨迹分析工具
	*/
	class GmxHandleFull
	{
	public:

		GmxHandleFull(int argc, char** argv);
		GmxHandleFull(const GmxHandleFull&) = delete;
		GmxHandleFull(GmxHandleFull&&) = delete;
		GmxHandleFull& operator=(const GmxHandleFull&) = delete;
		GmxHandleFull& operator=(GmxHandleFull&&) = delete;

		/**
		 * @brief 初始化，包括选择索引
		*/
		void init();

		/**
		 * @brief 读取第一帧数据
		 * @return 是否成功读取
		*/
		bool readFirstFrame();

		/**
		 * @brief 读取下一帧数据
		 * @return 是否成功读取
		*/
		bool readNextFrame();

		/**
		 * @brief 初始化一个空的可以储存分子位置的结构
		 * @param grp 选择的组号
		 * @return 空的结构
		*/
		void initPos(boxd& pos, int grp);

		/**
		 * @brief 初始化一个空的可以储存分子中心位置的结构
		 * @param grp 选择的组号
		 * @return 空的结构
		*/
		void initPosc(matd& posc, int grp);

		/**
		 * @brief 导入分子坐标信息
		 * @param pos 写入坐标的容器
		 * @param grp 组别
		*/
		void loadPosition(boxd& pos, int grp);

		/**
		 * @brief 导入分子坐标中心信息
		 * @param pos 写入坐标中心的容器
		 * @param grp 组别
		 * @param com 质心计算方式；质量(0)，几何(1)，电荷(2)
		*/
		void loadPositionCenter(matd& posc, int grp, int com = 0);

		/**
		 * @brief 得到一个文件句柄
		 * @param fnm 文件名
		 * @param writeInfo 是否写入额外文件信息，默认写入创建时间，创建路径等信息
		 * @return 文件句柄
		*/
		FILE* openWrite(std::string fnm, bool writeInfo = true);

		/**
		 * @brief 根据文件类型获取文件名
		 * @param ftp 文件类型, "efXVG", ...
		 * @return 文件名
		*/
		const char* get_ftp2fn(int ftp);

		/**
		 * @brief 根据文件类型获取文件名
		 * @param ftp 文件类型, "efXVG", ...
		 * @return 文件名
		*/
		const char* get_ftp2fn_null(int ftp);

		/**
		 * @brief 根据命令行选项获取文件名
		 * @param opt 命令行选项, "-o", ...
		 * @return 文件名
		*/
		const char* get_opt2fn(const char* opt);

		/**
		 * @brief 处理周期性问题
		 * @param dx dx
		 * @param box 盒子长度
		 * @return dx
		*/
		static double periodicity(double dx, double box);

	private:
		/**
		 * @brief 获取组中每个分子所含原子数
		 * @param grp 组别
		 * @return 每个分子中含原子数
		*/
		void get_natom_per_mol(int grp);

		/**
		 * @brief 获取组中分子每个原子的质量
		 * @param grp 组别
		 * @return 分子中每个原子质量数组
		*/
		void get_mass_charge(int grp);

	public:
		template <typename T>
		using stdMatrix = std::vector<std::vector<T>>;

		std::vector<const char*> desc = { "do some analysis for trajectort." };

		double               Lbox[3];      // Length of box 		
		char**               grpname;      // group names                
		int*                     ngx;      // sizes of groups             
		t_topology*              top;      // topology   
		t_inputrec*               ir;      // input rec, info from mdp file
		double                    dt;      // dt, ps
		double                  time;      // current time;
		double               preTime;      // pre time;
		int**                  index;      // indices for all groups     
		t_trxframe*               fr;      // frame information          
		int                    flags;      // flag . which (x, v, or f) to read, default(TRX_READ_X) 
		int                    ngrps;      // nr. of groups, default(1)               
		int                   nframe;      // nr. of frame                             
		std::vector<t_pargs>      pa;      // add some user-defined pargs             
		std::vector<t_filenm>    fnm;      // files for analysis

		stdMatrix<int>                   napm;        // number of atoms per molecule in each selection
		std::vector<int>                 nmol;        // number of nolecules in each selection 
		std::vector<stdMatrix<double>>   mass;        // mass of each group  
		std::vector<stdMatrix<double>>   charge;      // charge of each group 
		stdMatrix<double>                totMass;     // total mass of molecules in each group
		stdMatrix<double>                totCharge;   // total charge of molecules in each group
		
	protected:
		int                     argc;
		char**                  argv;
		t_trxstatus*          status;
		gmx_output_env_t*       oenv;
		int                     ePBC;
	}; // ! class GmxHandle

} // ! namespace itp
