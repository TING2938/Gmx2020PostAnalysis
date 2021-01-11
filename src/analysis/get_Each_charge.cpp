#include <itp/core>
#include <itp/fileio>

int main()
{
	auto d = itp::loadtxt("eachAtomCharge.dat");
	fmt::print("size of data is {} x {}.\n", d.rows(), d.cols());

	/* for each charge */
	Eigen::ArrayXXd dd = d.bottomRows(d.rows() / 2).colwise().mean();
	dd.resize(300, 36);

	Eigen::ArrayXXd CC_pos = dd.leftCols(18).rowwise().mean();
	Eigen::ArrayXXd CC_neg = dd.rightCols(18).rowwise().mean();

	std::ofstream file("AvgEachCharge.dat");
	for (int i = 0; i != CC_pos.size(); i++)
	{
		fmt::print(file, "{}\t{}\n", CC_pos(i), CC_neg(i));
	}

	/* for charge-time distribution */
	double atype[] = { 0.555621, -0.53234, 0.26482, 0.096227, 0.231484, -0.24232, 0.043224 };
	Eigen::ArrayXXd initCharge = itp::loadtxt("../init.dat");
	
	Eigen::ArrayXXd charge_pos(d.rows(), 7);
	Eigen::ArrayXXd charge_neg(d.rows(), 7);
	charge_pos.fill(0);
	charge_neg.fill(0);
	int count = 0;

	for (int i = 0; i != 7; ++i)
	{
		count = 0;
		for (int p = 0; p != 300; ++p)
		{
			if (std::abs(initCharge(p) - atype[i]) < 1e-3)
			{
				count++;

				for (int q = 0; q != 18; ++q)
				{
					charge_pos.col(i) += d.col(p + q * 300);
					charge_neg.col(i) += d.col(p + q * 300 + 5400);
				}
			}
		}
		
		charge_pos.col(i) /= (count * 18);
		charge_neg.col(i) /= (count * 18);
	}

	std::ofstream charge_time("charge_time.dat");
	Eigen::ArrayXXd charge_tot(d.rows(), 14);
	charge_tot << charge_pos, charge_neg;
	charge_time << charge_tot;


}
