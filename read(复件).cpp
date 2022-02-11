#include "atom.h"
#include "input.h"
#include "fun.h"
#include "atom.h"
#include "read.h"
#include <random>

Read::Read(){}
Read::~Read(){}
int Read::CountLines(string &filename, string A)
{
    ifstream ReadFile;
	string type;
	double x = 0.0;
	double y = 0.0;
	double z = 0.0;

    static int a[2];  
    int n = 0;
    int k = 0;
    string tmp; 
    ReadFile.open(filename,ios::in);
    if(ReadFile.fail())
    {
        return 0;
    }
    else
    {
        while(getline(ReadFile,tmp,'\n'))
        {
        	if(tmp.compare(A)==0)
		   {
			    k = n+1; 
			    a[0] = k; 
				//n = -1; 
			}
            n++;
            
        }
        ReadFile.close();
		}
        return a[0];
}

void Read::Init_position(Geo &geo)
 // read position from geo.in
 {
    ifstream ifs(INPUT.geo_dir);
	int linenumber1 = 0;
	int linenumber2 = 0;
	int PostionNumber = 0;	
	
	// initialize 
	
    linenumber1 = CountLines(INPUT.geo_dir,"%ATOMIC_POSTION");
	cout << linenumber1 << endl;
	//linenumber2 = CountLines(INPUT.geo_dir,"%ATOMIC_VELOCITY");
	//cout << linenumber2 << endl;
	//PostionNumber = linenumber2-linenumber1-1;
	for(int i = 1; i <= linenumber1;i++)
	{
		ifs.ignore(150, '\n');
	}
	for(int i = 0; i < INPUT.natom; i++)
	{
		ifs >> geo.atom[i].id >> geo.atom[i].x >> geo.atom[i].y >> geo.atom[i].z;
		periodicity_process(geo.atom[i].x, INPUT.cellpara1);
		periodicity_process(geo.atom[i].y, INPUT.cellpara2);
		periodicity_process(geo.atom[i].z, INPUT.cellpara3);
	}

	ifs.close();
 }

void Read::Init_volicity_fromfile(Geo &geo)
 // read position from geo.in
 {
    ifstream ifs(INPUT.geo_dir);
	int linenumber = 0; // Total lines of the file(geo.in)
	int linenumber1 = 0; // The number of rows in which "%ATOMIC_VELOCITY" resides
	int PostionNumber = 0;	// the number of lines we need read
	
	// initialize 
	linenumber1 = CountLines(INPUT.geo_dir,"%ATOMIC_VELOCITY");
	PostionNumber = INPUT.natom;	
	for(int i = 0; i < linenumber1;i++)
	{
		ifs.ignore(150, '\n');
	}
	for(int i = 0; i < PostionNumber; i++)
	{
		ifs >> geo.atom[i].id >> geo.atom[i].vx >> geo.atom[i].vy >> geo.atom[i].vz;
	}

	ifs.close();
 }

 void Read::Init_volicity_byRandomNumber(Geo &geo,double T)
 //Use random number to initialize velocity   T:temperature
 {
	const int SEED = 6666;
	default_random_engine random(SEED);
    std::uniform_real_distribution<double> dis(-0.5, 0.5);  
	for(int i = 0; i < INPUT.natom; i++)
	{
		geo.atom[i].vx = dis(random);
		geo.atom[i].vy = dis(random);
		geo.atom[i].vz = dis(random);
	}
	
	double m = 0.0; //The weight per atom m_i
	double v_cx = 0.0; //v_ix = v_ix - v_cx
	double v_cy = 0.0; //v_iy = v_iy - v_cy
	double v_cz = 0.0; //v_iz = v_iz - v_cz
	double hx1 = 0.0; // h1 = sum of m_i*v_ix
	double hy1 = 0.0;
	double hz1 = 0.0;  
	//double h2 = 0.0; // h2 = sum of m_i
	double factor = 0.0;
	double factory = 0.0;
	double factorz = 0.0;
	double fx = 0.0; // sum of (v'_ix)^2
	double fy = 0.0; 
	double fz = 0.0; 
	double K = 1.380649; //e-23Boltzmann constant
	double NA = 6.02214076;//e23
	double EV  = 1.602176634; //e-19
	//m = INPUT.mass*1e-3/6.02214086e23;
	m = INPUT.mass*1e-1/NA; 
	//h2 = m*INPUT.natom;
	/*
	for(int i = 0; i < INPUT.natom; i++)
	{
		hx1 = hx1+m*geo.atom[i].vx;
		hy1 = hy1+m*geo.atom[i].vy;
		hz1 = hz1+m*geo.atom[i].vz;
	}
	v_cx = hx1/h2;
	v_cy = hy1/h2;
	v_cz = hz1/h2;
	*/
	for(int i = 0; i < INPUT.natom; i++)
	{
		hx1 = hx1+geo.atom[i].vx;
		hy1 = hy1+geo.atom[i].vy;
		hz1 = hz1+geo.atom[i].vz;
	}	
	v_cx = hx1/INPUT.natom;
	v_cy = hy1/INPUT.natom;
	v_cz = hz1/INPUT.natom;
	for(int i = 0; i < INPUT.natom; i++)
	{
		geo.atom[i].vx = geo.atom[i].vx - v_cx;
		geo.atom[i].vy = geo.atom[i].vy - v_cy;
		geo.atom[i].vz = geo.atom[i].vz - v_cz;
		/*
		cout << i << endl;
		cout << geo.atom[i].vx << endl; 
		cout << geo.atom[i].vy << endl;
		cout << geo.atom[i].vz << endl;
		*/
	}

	for(int i = 0; i < INPUT.natom; i++)
	{
		fx = fx+pow(geo.atom[i].vx,2);
		fy = fy+pow(geo.atom[i].vy,2);
		fz = fz+pow(geo.atom[i].vz,2);
	}
	//cout << fx <<" " << fy <<" " << fz <<endl;
	//cout << m << " "<< K <<" "<< INPUT.natom<<" "<< T <<endl;
	factor = pow(3.0*INPUT.natom*K*T/((INPUT.mass / NA *1e-3)*(fx+fy+fz)*1e4),0.5);

	
	for(int i = 0; i < INPUT.natom; i++)
	{
		geo.atom[i].vx = factor*geo.atom[i].vx;
		geo.atom[i].vy = factor*geo.atom[i].vy;
		geo.atom[i].vz = factor*geo.atom[i].vz;
	}
	
	return;	
 }
