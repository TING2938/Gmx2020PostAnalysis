[](#起 "起")起[](#起)
-----------------

[Gromacs](http://www.gromacs.org/)是一款十分优秀的分子动力学软件，其并行化运行速度快，且本身自带丰富的后处理命令。但在科研过程中难免会出现自带分析工具难以满足我们需求的情况，这时候就需要我们自己针对Gromacs产生的轨迹文件编写符合我们需求的后处理程序。本文介绍了针对Gromacs-2020系列版本的轨迹后处理分析过程，程序编写基于本人在研究生期间开发的后处理框架，开发平台为Windows 10-Visual Studio Community 2019，在Windows上编写代码就可以在Windows上编译运行，也可以把代码Copy到Linux服务器上编译运行。

### [](#准备文件和工具 "准备文件和工具")准备文件和工具[](#准备文件和工具)

1.  Gmx2020后处理框架
    
2.  [Visual Studio Community 2019](https://visualstudio.microsoft.com/zh-hans/vs/)
    
    安装时只需勾选`使用c++的桌面开发`
    

准备完成后，解压缩后处理框架，进入`Gmx2020_PostAnalysis_for_VS2019/Gmx2020PostAnalysis/`目录中，双击`Gmx2020PostAnalysis.sln`即可打开VS2019工程。

[](#承 "承")承[](#承)
-----------------

### [](#简单程序代码示例 "简单程序代码示例")简单程序代码示例[](#简单程序代码示例)

后处理程序由一个简单的示例引入：

```c++
/**
 * @file template.cpp
 * @author yeting (yeting2938@hust.edu.cn)
 * @brief a simple demo of gromacs trajectry analysis tool
 * @version 0.1
 * @date 2019/05/05
 * 
 * @copyright Copyright (c) 2019
 * 
 */


// 头文件，包含常用库和函数
#include <itp/gmx>

/**
 * @brief 位于命名空间`itp`下的类`GmxHandle`里面定义了绝大多数在轨迹分析过程中
 * 会用到的数据以及函数，比如体系的拓扑结构、分子与原子的数目、分子的质量与电量
 * 等数据结构，获取分子位置、速度、受力等函数，会在下节详细讲解。这里是构造一个
 * 继承类，里面可以定义一些在本代码中会用到的数据以及函数。
 */
class Handle : public itp::GmxHandle
{
public:
	using GmxHandle::GmxHandle;

public:
	itp::boxd pos1;
	itp::matd posc1;
};


// `gmx_main`是一个宏，其参数`temp`是一个函数名称，可以把`gmx_main`当作`main`函数。
gmx_main(temp)
{
	// 实例化一个`Handle`对象，其参数`argc`与`argv`定义在宏`gmx_main`里面
	Handle hd(argc, argv);

	/**
	 * @brief 该后处理框架对体系的分析基于`Gromacs`中组(Group)的概念，即仅对
	 * 用户选择的组进行分析，选择的组里面包含了该组中所有原子的序号。这里可以定
	 * 义该程序所需组的个数，然后为每一个组定义一个变量方便后面代码操作。
	 */
	hd.ngrps = 1; // number of group(s);
	int grp = 0;  // selected group;

	/**
	 * @brief add some user-defined pargs. here
	 * 在这里定义由命令行输入的参数，采用`Gromacs`自带参数解析系统，支持的变量类型有：
	 *         `int, int64, real, time, str, bool, rvec, enum` 
	 * 详情可在`Gromacs`源代码文件`include/gromacs/commandline/pargs.h`中查看
	 * @note 浮点数要用`real`类型，不能直接用`double`类型
	 */
	int nbin = 100;
	int dim = 2;     // 0(x), 1(y), 2(z);
	real lowPos = 0; // nm;
	real upPos = 30; // nm;

	// 在这里编写命令行参数接口
	hd.pa = {
		{ "-nbin", FALSE, etINT,  {&nbin}, "nbins."},
		{ "-dim", FALSE, etINT, {&dim}, "dim, 0(x), 1(y), 2(z)"},
		{ "-up", FALSE, etREAL, {&upPos}, "up position of region of molecule/ion (nm)" },
		{ "-low", FALSE, etREAL, {&lowPos}, "low position of region of molecule/ion (nm)" }
	};
	
	// 在这里定义后处理过程中涉及到的输入/输出文件名称
	hd.fnm = {
		{ efXVG, "-o", "analysisOutput", ffWRITE }
	};

	/**
	 * @brief 初始化，在这个函数内会处理命令行选项，选择分析组，分析拓扑信息，并
	 * 计算组内原子的电荷量与质量。
	 * @note 该函数有一个可选参数，默认为`true`，代表选择的组内具有同种分子，否则为`false`。
	 */
	hd.init();

	// 读取第一帧轨迹数据。数据更新在`hd.fr`内。
	hd.readFirstFrame();

	/**
	 * @brief `hd.initPos(grp)`与`hd.initPosc(grp)`两个函数初始化两个矩阵类型数据
	 * 结构，分别储存分子中原子轨迹信息与分子中心轨迹信息。函数的参数为选择组的组号。
	 * 	`hd.initPos(grp)`返回的数据大小为: [分子个数]x[每个分子中原子个数]x[3]
	 * 	`hd.initPosc(grp)`返回的数据大小为: [分子个数]x[3]
	 */
	hd.pos1 = hd.initPos(grp);
	hd.posc1 = hd.initPosc(grp);

	/**
	 * @brief `itp::vecd`为一种数组类型，数组中每个元素类型为`double`，数组
	 * 长度为`nbin`。函数`fill`指向数组中填充相同的数值。
	 * 相似的数据结构还有：
	 * 		`itp::veci`: `int`型数组
	 * 		`itp::vecd`: `double`型数组
	 * 		`itp::mati`: `int`型矩阵
	 * 		`itp::matd`: `double`型矩阵
	 * 		`itp::boxd`: `double`型三维矩阵
	 * 这些类型为C++开源矩阵库Eigen(https://eigen.tuxfamily.org)中类型的别名
	 */
	itp::vecd density(nbin);
	density.fill(0);

	double dbin = (upPos - lowPos) / nbin;
	int me = 0;

	/* ------------- main loop for each frame ----------------- */
	do
	{
		/**
		 * @brief 获取位置以及位置中心，并储存在相应的数据结构中，
		 * `loadPositionCenter`还有第三个可选参数，默认为0（代表取质量中心）
		 * 还有1（代表几何中心）和2（代表电荷中心）
		 */
		hd.loadPosition(hd.pos1, grp);
		hd.loadPositionCenter(hd.posc1, grp);

		for (int i = 0; i < hd.nmol[grp]; i++)
		{
			if (lowPos <= hd.posc1(i, dim) && hd.posc1(i, dim) < upPos)
			{
				me = (hd.posc1(i, dim) - lowPos) / dbin;
				density[me]++;
			}
		}
	} while (hd.readNextFrame());
	/* --------------- main loop end -------------------- */

	density /= hd.nframe;
	
	/* ------------- for output --------------- */

	/**
	 * @brief 打开一个文件，文件名由命令行参数"-o"决定
	 * `openWrite`函数还有第二个可选参数，默认为`true`，代表在文件开始处写入
	 * 文件生成的相关信息（如路径、创建时间、命令行参数等）
	 */
	auto ofile = hd.openWrite(hd.get_opt2fn("-o"));

	for (int i = 0; i < nbin; i++)
	{
		/**
		 * @brief 输出数据到文件。这里输出操作采用的是一个
		 * 开源库{fmt}(https://github.com/fmtlib/fmt)
		 */
		fmt::print(ofile, "{:8.5f} {10.5f}\n", lowPos + dbin / 2 + dbin*i, density[i]);
	}
	fclose(ofile);

	return 0;
}
```

以上程序统计指定范围内（`-low`与`-up`之间）、指定维度（`-dim`）分子的个数分布。

[](#转 "转")转[](#转)
-----------------

*   class `itp::GmxHandle`
    
    位于命名空间`itp`下的类`GmxHandle`里面定义了绝大多数在轨迹分析过程中会用到的数据以及函数，比如体系的拓扑结构、分子与原子的数目、分子的质量与电量等数据结构，获取分子位置、速度、受力等函数
    
    *   数据成员
        
        *   **desc** (std::vector<const char\*>)
            
            描述本程序用途的字符串
            
        *   **Lbox** (double\[3\])
            
            模拟盒子的长度，每当采用loadPosition等函数读取轨迹都会更新
            
        *   **grpname** (char\*\*)
            
            选择组的名称
            
        *   **ngx** (int\*)
            
            选择组的中含有多少原子
            
        *   **top** (t\_topology\*)
            
            拓扑信息。里面包含所有拓扑信息
            
        *   **ir** (t\_inputrec\*)
            
            输入参数，包含mdp文件中所有信息
            
        *   **dt** (double)
            
            模拟时间步长，ps
            
        *   **time** (double)
            
            模拟的时刻
            
        *   **preTime** (double)
            
            模拟的上一步时刻
            
        *   **index** (int\*\*)
            
            原子索引，包含选择的每一组中的原子序号
            
        *   **fr** (t\_trxframe\*)
            
            轨迹相关信息
            
        *   **flags** (int)
            
            flags，控制是否读取位置、速度或力，默认为（TRX\_READ\_X）
            
        *   **ngrps** (int)
            
            选择的组的数目，默认为1
            
        *   **nframe** (int)
            
            读取到的轨迹总帧数
            
        *   **pa** (std::vector)
            
            储存命令行参数
            
        *   **fnm** (std::vector)
            
            储存命令行文件输入/输出
            
        *   **napm** (veci)
            
            在选择的组中分子中包含的原子数
            
        *   **nmol** (veci)
            
            在选择的组中的分子数
            
        *   **mass** (Vec)
            
            在选择的组中一个分子中每个原子的质量
            
        *   **charge** (Vec)
            
            在选择的组中一个分子中每个原子的电荷量
            
    *   成员函数
        
        *   `GmxHandle(int argc, char** argv)`
            
            构造函数
            
            参数：`argc`与`argv`为`main()`函数的参数
            
        *   `void init(bool bFullMolecule=true)`
            
            初始化，在这个函数内会处理命令行选项，选择分析组，分析拓扑信息，并计算组内原子的电荷量与质量。
            
            注：可选参数`bFullMolecule`，默认为`true`，代表选择的组内具有同种分子，否则为`false`。
            
        *   `bool readFirstFrame()`
            
            读取第一帧轨迹数据，
            
            返回：是否读取成功
            
        *   `bool readNextFrame()`
            
            读取下一帧轨迹数据
            
            返回：是否读取成功
            
        *   `boxd initPos(int grp)`
            
            初始化一个空的可以储存分子位置信息的数据结构
            
            参数：`grp`为组号
            
        *   `matd initPosc(int grp)`
            
            初始化一个空的可以储存分子中心位置的结构
            
            参数：`grp`为组号
            
        *   `void loadPosition(boxd& pos, int grp)`
            
            读取位置信息，并储存到结构`pos`中，`grp`为组号
            
        *   `void loadPositionCenter(matd& posc, int grp, int com=0)`
            
            读取位置中心信息，并储存到结构`pos`中，`grp`为组号
            
        *   `void loadVelocity(boxd& vel, int grp)`
            
            读取速度信息，并储存到结构`vel`中，`grp`为组号
            
        *   `void loadVelocityCenter(matd& velc, int grp)`
            
            读取速度中心信息，并储存到结构`velc`中，`grp`为组号
            
        *   `void loadForce(boxd& force, int grp)`
            
            读取力信息，并储存到结构`force`中，`grp`为组号
            
        *   `void loadForceCenter(matd& forcec, int grp)`
            
            读取力中心信息，并储存到结构`forcec`中，`grp`为组号
            
        *   `void loadLJParameter(int grp1, int grp2, matd& c6, matd& c12)`
            
            读取LJ参数信息，储存到`c6`和`c12`中，两组组号分别为`grp1`和`grp2`
            
        *   `FILE* openWrite(std::string fnm, bool writeInfo=true)`
            
            讲注：打开一个文件
            
            参数：
            
            `fnm`：文件名
            
            `writeInfo`： 可选参数，默认为`true`，代表在文件开始处写入文件生成的相关信息（如路径、创建时间、命令行参数等）
            
        *   `const char* get_ftp2fn(int ftp)`
            
            通过命令行参数类型获取文件名称
            
            参数：
            
            `ftp`：命令行参数类型，有`efXVG`、`efNDX`、`efTPR`、`efNDX`等值可选，可在文件`include/gromacs/fileio/filetypes.h`中查看完整可选范围。
            
        *   `const char* get_ftp2fn_null(int ftp)`
            
            同上
            
        *   `const char* get_opt2fn(const char* opt)`
            
            通过命令行选项获取文件名称
            
            参数：
            
            `opt`：命令行文件输入/输出选项
            
        *   `static double periodicity(double dx, double box)`
            
            静态成员函数，处理周期性问题
            
            参数：
            
            `dx`：两位置之差
            
            `box`：在计算`dx`维度上的盒子长度
            

*   此后处理框架中用到的外部开源库
    
    *   [Eigen](https://eigen.tuxfamily.org/)
        
        一个开源的线性代数库，里面有方便使用的矩阵、向量等类型用法参照[Eigen与Matlab操作对照](https://ting2938.github.io/2020/07/18/Eigen%E4%B8%8EMatlab%E6%93%8D%E4%BD%9C%E5%AF%B9%E7%85%A7)
        
        在本框架中常用的数据结构有：
        
        `itp::veci`: `int`型数组
        
        `itp::vecd`: `double`型数组
        
        `itp::mati`: `int`型矩阵
        
        `itp::matd`: `double`型矩阵
        
        `itp::boxd`: `double`型三维矩阵
        
    *   [fmt](https://github.com/fmtlib/fmt)
        
        一个开源的C++输出格式控制库，用法类似`Python`的格式化字符串方式，如
        
        ```null
        
        fmt::print("Hello {}!\n", "world");
        fmt::print("4.2 + 2.1 = {:5.2f}", 6.3);
        
        
        FILE* ofile("output.dat", "w");
        fmt::print(ofile, "Hello {}!\n", "world");
        std::ofstream ofsfile("output.dat");
        fmt::print(ofsfile, "Hello {}!\n", "world");
        ```
        

[](#合 "合")合[](#合)
-----------------

运用本框架进行Gromacs后处理分析过程中遇到问题可以联系我

> Name : 叶挺
> 
> Email: [yeting2938@hust.edu.cn](mailto:yeting2938@hust.edu.cn)

* * *

## How to compile

```bash
export CODE_INSTALL_PATH=~/software/myCode  # where to install 
./INSTALL src/general/calc_energy.cpp # path to your source code
```
