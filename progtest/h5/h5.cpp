#include  "Alat/armadillo.hpp"

void check(std::string file1, std::string file2)
{
	std::string cmd = "h5dump "+ file1 + ".h5 >> " + file1;
	system(cmd.c_str());	
	cmd = "h5dump "+ file2 + ".h5 >> " + file2;
	system(cmd.c_str());	
	cmd = "diff  -y -W 70 "+ file1 + " " + file2;
	system(cmd.c_str());	
}

int main()
{
	arma::mat v(3,2,arma::fill::randu), w;
	std::string file1("file1"), file2("file2");
	arma::hdf5_name spec1(file1+".h5","TOTO/TATA");
	arma::hdf5_name spec2(file2+".h5","TOTO/TATA");
	v.save(spec1);
	w.load(spec1);
	w.save(spec2);
	check(file1, file2);
	return 0;
}