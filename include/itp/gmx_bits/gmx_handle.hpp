
namespace itp
{
	/**
	 * @brief Gromacs轨迹分析工具
	*/
	class GmxHandle
	{
	public:

		GmxHandle(int argc, char** argv);
		GmxHandle(const GmxHandle&) = delete;
		GmxHandle(GmxHandle&&) = delete;
		GmxHandle& operator=(const GmxHandle&) = delete;
		GmxHandle& operator=(GmxHandle&&) = delete;

		/**
		 * @brief 初始化，包括选择索引
		*/
		void init(bool fullMolecule = true);

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
		boxd initPos(int grp);

		/**
		 * @brief 初始化一个空的可以储存分子中心位置的结构
		 * @param grp 选择的组号
		 * @return 空的结构
		*/
		matd initPosc(int grp);

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
		 * @brief 导入分子速度信息
		 * @param pos 写入速度的容器
		 * @param grp 组别
		*/
		void loadVelocity(boxd& vel, int grp);

		/**
		 * @brief 导入分子速度中心信息
		 * @param pos 写入速度中心的容器
		 * @param grp 组别
		 * @param com 质心计算方式；质量(0)，几何(1)，电荷(2)
		*/
		void loadVelocityCenter(matd& velc, int grp, int com = 0);

		/**
		 * @brief 导入分子速度信息
		 * @param pos 写入速度的容器
		 * @param grp 组别
		*/
		void loadForce(boxd& force, int grp);

		/**
		 * @brief 导入分子速度中心信息
		 * @param pos 写入速度中心的容器
		 * @param grp 组别
		 * @param com 质心计算方式；质量(0)，几何(1)，电荷(2)
		*/
		void loadForceCenter(matd& forcec, int grp, int com = 0);

		/**
		 * @brief 获取两组分子之间原子的LJ参数
		 * @param group1 组别1
		 * @param group2 组别2
		 * @param c6 c6参数
		 * @param c12 c12参数
		*/
		void loadLJParameter(int group1, int group2, matd& c6, matd& c12);

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
		int get_natom_per_mol(int grp);

		/**
		 * @brief 获取组中分子每个原子的质量
		 * @param grp 组别
		 * @return 分子中每个原子质量数组
		*/
		vecd get_mass(int grp);

		/**
		 * @brief 获取组中分子每个原子的电荷量
		 * @param grp 组别
		 * @return 分子中每个原子电荷量数组
		*/
		vecd get_charge(int grp);

	public:
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

		veci                    napm;      // number of atoms per molecule in each selection
		veci                    nmol;      // number of nolecules in each selection 
		Vec<vecd>               mass;      // mass of each group  
		Vec<vecd>             charge;      // charge of each group 
		
	protected:
		int                     argc;
		char** argv;
		t_trxstatus* status;
		gmx_output_env_t* oenv;
		int                     ePBC;
	}; // ! class GmxHandle

} // ! namespace itp
