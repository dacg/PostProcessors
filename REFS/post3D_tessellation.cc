//################################################
//############ Tessellation generator ############
//########## Based in Voro++ Examples ############
//############## D. Cantor 16/10/14 ##############
//################################################
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <iomanip>
#include <cmath>
#include <time.h>
#include "voro++.hh"

using namespace std;
using namespace voro;

// Setting a function that returns a random double between 0 and 1

double rnd() {return double(rand())/RAND_MAX;}


// Setting constants
// This is the number of blocks in which the container is divided.
const int n_x=8,n_y=8,n_z=8;

int main() {
	
	srand (time(NULL));
	
	// Reading the main input file
	
	string line;
	ifstream file("./POSTPRO/VORONOI/0main");
	
	int n_particles;
	int n_sub_par;
	double a_voro;
	
	if (file.is_open())
	{	
		// Reading the number of particles
		getline(file,line);
		string s_n_particles = line.substr(0, line.find("#", 0)); 
		
		n_particles = stoi(s_n_particles);
		
		// Reading the number of sub particles 
		getline(file,line);
		string s_n_sub_par = line.substr(0, line.find("#", 0)); 
		
		n_sub_par = stoi(s_n_sub_par);
		
		// Reading the randomness parameter
		getline(file,line);
		string s_a_voro = line.substr(0, line.find("#", 0)); 
		
		a_voro = stod(s_a_voro);
		
	}
	file.close();
	
	// Variables for the lecture of the faces of the particle
	int parsing = 0;
	string input_file;
	
 	double x_min;
 	double y_min;
	double z_min;
	double x_max;
 	double y_max;
	double z_max;
	
	cout << "########################################" << "\n";
	cout << "      Setting up the tessellation       " << "\n";
	
	// For each particle
	while (parsing < n_particles)
	{
		parsing++;
		input_file = "./POSTPRO/VORONOI/input" + to_string(parsing);
		ifstream file_imp(input_file);
		
		// Reading the maximum and minimum vertices of the particle
		string temp;
		
		//Reading the minimum values
		getline(file_imp,temp);
		x_min = stod(temp);
		
		getline(file_imp,temp);
		y_min = stod(temp);
		
		getline(file_imp,temp);
		z_min = stod(temp);
		
		// Reading the maximum values
		getline(file_imp,temp);
		x_max = stod(temp);
		
		getline(file_imp,temp);
		y_max = stod(temp);
		
		getline(file_imp,temp);
		z_max = stod(temp);
		
		// Reading the number of faces
		getline(file_imp,temp);
		int n_faces = stoi(temp);
		
		double array[n_faces][4];
		
		// For each face
		for (int count = 0; count<n_faces; count++)
		{
			// Reading the data and splitting it up 
			getline(file_imp,temp);
			string s_dist = temp.substr(0,16);
			string s_x = temp.substr(19,35);
			string s_y = temp.substr(38,54);
			string s_z = temp.substr(57,73); 
			
			// Storing the data of walls
			array[count][0] = stod(s_dist);
			array[count][1] = stod(s_x);
			array[count][2] = stod(s_y);
			array[count][3] = stod(s_z);
			
		}
		
		// Creating a container with the geometry given above, and make it
		// non-periodic in each of the three coordinates. Allocate space for 8
		// particles within each computational block.
		container con(x_min,x_max,y_min,y_max,z_min,z_max,n_x,n_y,n_z,
		false,false,false,8);
		
		// Declaring walls
		// This step is extremely inefficient and should be improved. The origin of this
		// inefficiency is the imposibility to name variables dynamically in c++.
		// Attention!! The following code is done for a maximum of 130 faces.
		// If the number of faces of the polyhedron exceeds this number the program will
		// stop and will generate a commentary in terminal.
		
		if (n_faces>130)
		{
			cout << "##################################################" << "\n";
			cout << "##################################################" << "\n";
			cout << "############# TESSELATION FAILED #################" << "\n";
			cout << "########## N_FACES > 130. See source #############" << "\n";
			cout << "##################################################" << "\n";
			cout << "##################################################" << "\n";
			return 0;
		}
		
		// Declaring and adding the corresponding walls for each face
		for (int count = 0; count<n_faces; count++)
		{	
			if(count==0)
			{
				wall_plane p0(array[count][1],array[count][2],array[count][3],array[count][0]);
				con.add_wall(p0);
				continue;
			}
			if(count==1)
			{
				wall_plane p1(array[count][1],array[count][2],array[count][3],array[count][0]);
				con.add_wall(p1);
				continue;
			}
			if(count==2)
			{
				wall_plane p2(array[count][1],array[count][2],array[count][3],array[count][0]);
				con.add_wall(p2);
				continue;
			}
			if(count==3)
			{
				wall_plane p3(array[count][1],array[count][2],array[count][3],array[count][0]);
				con.add_wall(p3);
				continue;
			}
			if(count==4)
			{
				wall_plane p4(array[count][1],array[count][2],array[count][3],array[count][0]);
				con.add_wall(p4);
				continue;
			}
			if(count==5)
			{
				wall_plane p5(array[count][1],array[count][2],array[count][3],array[count][0]);
				con.add_wall(p5);
				continue;
			}
			if(count==6)
			{
				wall_plane p6(array[count][1],array[count][2],array[count][3],array[count][0]);
				con.add_wall(p6);
				continue;
			}
			if(count==7)
			{
				wall_plane p7(array[count][1],array[count][2],array[count][3],array[count][0]);
				con.add_wall(p7);
				continue;
			}
			if(count==8)
			{
				wall_plane p8(array[count][1],array[count][2],array[count][3],array[count][0]);
				con.add_wall(p8);
				continue;
			}
			if(count==9)
			{
				wall_plane p9(array[count][1],array[count][2],array[count][3],array[count][0]);
				con.add_wall(p9);
				continue;
			}
			if(count==10)
			{
				wall_plane p10(array[count][1],array[count][2],array[count][3],array[count][0]);
				con.add_wall(p10);
				continue;
			}
			if(count==11)
			{
				wall_plane p11(array[count][1],array[count][2],array[count][3],array[count][0]);
				con.add_wall(p11);
				continue;
			}
			if(count==12)
			{
				wall_plane p12(array[count][1],array[count][2],array[count][3],array[count][0]);
				con.add_wall(p12);
				continue;
			}
			if(count==13)
			{
				wall_plane p13(array[count][1],array[count][2],array[count][3],array[count][0]);
				con.add_wall(p13);
				continue;
			}
			if(count==14)
			{
				wall_plane p14(array[count][1],array[count][2],array[count][3],array[count][0]);
				con.add_wall(p14);
				continue;
			}
			if(count==15)
			{
				wall_plane p15(array[count][1],array[count][2],array[count][3],array[count][0]);
				con.add_wall(p15);
				continue;
			}
			if(count==16)
			{
				wall_plane p16(array[count][1],array[count][2],array[count][3],array[count][0]);
				con.add_wall(p16);
				continue;
			}
			if(count==17)
			{
				wall_plane p17(array[count][1],array[count][2],array[count][3],array[count][0]);
				con.add_wall(p17);
				continue;
			}
			if(count==18)
			{
				wall_plane p18(array[count][1],array[count][2],array[count][3],array[count][0]);
				con.add_wall(p18);
				continue;
			}
			if(count==19)
			{
				wall_plane p19(array[count][1],array[count][2],array[count][3],array[count][0]);
				con.add_wall(p19);
				continue;
			}
			if(count==20)
			{
				wall_plane p20(array[count][1],array[count][2],array[count][3],array[count][0]);
				con.add_wall(p20);
				continue;
			}
			if(count==21)
			{
				wall_plane p21(array[count][1],array[count][2],array[count][3],array[count][0]);
				con.add_wall(p21);
				continue;
			}
			if(count==22)
			{
				wall_plane p22(array[count][1],array[count][2],array[count][3],array[count][0]);
				con.add_wall(p22);
				continue;
			}
			if(count==23)
			{
				wall_plane p23(array[count][1],array[count][2],array[count][3],array[count][0]);
				con.add_wall(p23);
				continue;
			}
			if(count==1)
			{
				wall_plane p1(array[count][1],array[count][2],array[count][3],array[count][0]);
				con.add_wall(p1);
				continue;
			}
			if(count==24)
			{
				wall_plane p24(array[count][1],array[count][2],array[count][3],array[count][0]);
				con.add_wall(p24);
				continue;
			}
			if(count==25)
			{
				wall_plane p25(array[count][1],array[count][2],array[count][3],array[count][0]);
				con.add_wall(p25);
				continue;
			}
			if(count==26)
			{
				wall_plane p26(array[count][1],array[count][2],array[count][3],array[count][0]);
				con.add_wall(p26);
				continue;
			}
			if(count==27)
			{
				wall_plane p27(array[count][1],array[count][2],array[count][3],array[count][0]);
				con.add_wall(p27);
				continue;
			}
			if(count==28)
			{
				wall_plane p28(array[count][1],array[count][2],array[count][3],array[count][0]);
				con.add_wall(p28);
				continue;
			}
			if(count==29)
			{
				wall_plane p29(array[count][1],array[count][2],array[count][3],array[count][0]);
				con.add_wall(p29);
				continue;
			}
			if(count==30)
			{
				wall_plane p30(array[count][1],array[count][2],array[count][3],array[count][0]);
				con.add_wall(p30);
				continue;
			}
			if(count==31)
			{
				wall_plane p31(array[count][1],array[count][2],array[count][3],array[count][0]);
				con.add_wall(p31);
				continue;
			}
			if(count==32)
			{
				wall_plane p32(array[count][1],array[count][2],array[count][3],array[count][0]);
				con.add_wall(p32);
				continue;
			}
			if(count==33)
			{
				wall_plane p33(array[count][1],array[count][2],array[count][3],array[count][0]);
				con.add_wall(p33);
				continue;
			}
			if(count==34)
			{
				wall_plane p34(array[count][1],array[count][2],array[count][3],array[count][0]);
				con.add_wall(p34);
				continue;
			}
			if(count==35)
			{
				wall_plane p35(array[count][1],array[count][2],array[count][3],array[count][0]);
				con.add_wall(p35);
				continue;
			}
			if(count==36)
			{
				wall_plane p36(array[count][1],array[count][2],array[count][3],array[count][0]);
				con.add_wall(p36);
				continue;
			}
			if(count==37)
			{
				wall_plane p37(array[count][1],array[count][2],array[count][3],array[count][0]);
				con.add_wall(p37);
				continue;
			}
			if(count==38)
			{
				wall_plane p38(array[count][1],array[count][2],array[count][3],array[count][0]);
				con.add_wall(p38);
				continue;
			}
			if(count==39)
			{
				wall_plane p39(array[count][1],array[count][2],array[count][3],array[count][0]);
				con.add_wall(p39);
				continue;
			}
			if(count==40)
			{
				wall_plane p40(array[count][1],array[count][2],array[count][3],array[count][0]);
				con.add_wall(p40);
				continue;
			}
			if(count==41)
			{
				wall_plane p41(array[count][1],array[count][2],array[count][3],array[count][0]);
				con.add_wall(p41);
				continue;
			}
			if(count==42)
			{
				wall_plane p42(array[count][1],array[count][2],array[count][3],array[count][0]);
				con.add_wall(p42);
				continue;
			}
			if(count==43)
			{
				wall_plane p43(array[count][1],array[count][2],array[count][3],array[count][0]);
				con.add_wall(p43);
				continue;
			}
			if(count==44)
			{
				wall_plane p44(array[count][1],array[count][2],array[count][3],array[count][0]);
				con.add_wall(p44);
				continue;
			}
			if(count==45)
			{
				wall_plane p45(array[count][1],array[count][2],array[count][3],array[count][0]);
				con.add_wall(p45);
				continue;
			}
			if(count==46)
			{
				wall_plane p46(array[count][1],array[count][2],array[count][3],array[count][0]);
				con.add_wall(p46);
				continue;
			}
			if(count==47)
			{
				wall_plane p47(array[count][1],array[count][2],array[count][3],array[count][0]);
				con.add_wall(p47);
				continue;
			}
			if(count==48)
			{
				wall_plane p48(array[count][1],array[count][2],array[count][3],array[count][0]);
				con.add_wall(p48);
				continue;
			}
			if(count==49)
			{
				wall_plane p49(array[count][1],array[count][2],array[count][3],array[count][0]);
				con.add_wall(p49);
				continue;
			}
			if(count==50)
			{
				wall_plane p50(array[count][1],array[count][2],array[count][3],array[count][0]);
				con.add_wall(p50);
				continue;
			}
			if(count==51)
			{
				wall_plane p51(array[count][1],array[count][2],array[count][3],array[count][0]);
				con.add_wall(p51);
				continue;
			}
			if(count==52)
			{
				wall_plane p52(array[count][1],array[count][2],array[count][3],array[count][0]);
				con.add_wall(p52);
				continue;
			}
			if(count==53)
			{
				wall_plane p53(array[count][1],array[count][2],array[count][3],array[count][0]);
				con.add_wall(p53);
				continue;
			}
			if(count==54)
			{
				wall_plane p54(array[count][1],array[count][2],array[count][3],array[count][0]);
				con.add_wall(p54);
				continue;
			}
			if(count==55)
			{
				wall_plane p55(array[count][1],array[count][2],array[count][3],array[count][0]);
				con.add_wall(p55);
				continue;
			}
			if(count==56)
			{
				wall_plane p56(array[count][1],array[count][2],array[count][3],array[count][0]);
				con.add_wall(p56);
				continue;
			}
			if(count==57)
			{
				wall_plane p57(array[count][1],array[count][2],array[count][3],array[count][0]);
				con.add_wall(p57);
				continue;
			}
			if(count==58)
			{
				wall_plane p58(array[count][1],array[count][2],array[count][3],array[count][0]);
				con.add_wall(p58);
				continue;
			}
			if(count==59)
			{
				wall_plane p59(array[count][1],array[count][2],array[count][3],array[count][0]);
				con.add_wall(p59);
				continue;
			}
			if(count==60)
			{
				wall_plane p60(array[count][1],array[count][2],array[count][3],array[count][0]);
				con.add_wall(p60);
				continue;
			}
			if(count==61)
			{
				wall_plane p61(array[count][1],array[count][2],array[count][3],array[count][0]);
				con.add_wall(p61);
				continue;
			}
			if(count==62)
			{
				wall_plane p62(array[count][1],array[count][2],array[count][3],array[count][0]);
				con.add_wall(p62);
				continue;
			}
			if(count==63)
			{
				wall_plane p63(array[count][1],array[count][2],array[count][3],array[count][0]);
				con.add_wall(p63);
				continue;
			}
			if(count==64)
			{
				wall_plane p64(array[count][1],array[count][2],array[count][3],array[count][0]);
				con.add_wall(p64);
				continue;
			}
			if(count==65)
			{
				wall_plane p65(array[count][1],array[count][2],array[count][3],array[count][0]);
				con.add_wall(p65);
				continue;
			}
			if(count==66)
			{
				wall_plane p66(array[count][1],array[count][2],array[count][3],array[count][0]);
				con.add_wall(p66);
				continue;
			}
			if(count==67)
			{
				wall_plane p67(array[count][1],array[count][2],array[count][3],array[count][0]);
				con.add_wall(p67);
				continue;
			}
			if(count==68)
			{
				wall_plane p68(array[count][1],array[count][2],array[count][3],array[count][0]);
				con.add_wall(p68);
				continue;
			}
			if(count==69)
			{
				wall_plane p69(array[count][1],array[count][2],array[count][3],array[count][0]);
				con.add_wall(p69);
				continue;
			}
			if(count==70)
			{
				wall_plane p70(array[count][1],array[count][2],array[count][3],array[count][0]);
				con.add_wall(p70);
				continue;
			}
			if(count==71)
			{
				wall_plane p71(array[count][1],array[count][2],array[count][3],array[count][0]);
				con.add_wall(p71);
				continue;
			}
			if(count==72)
			{
				wall_plane p72(array[count][1],array[count][2],array[count][3],array[count][0]);
				con.add_wall(p72);
				continue;
			}
			if(count==73)
			{
				wall_plane p73(array[count][1],array[count][2],array[count][3],array[count][0]);
				con.add_wall(p73);
				continue;
			}
			if(count==74)
			{
				wall_plane p74(array[count][1],array[count][2],array[count][3],array[count][0]);
				con.add_wall(p74);
				continue;
			}
			if(count==75)
			{
				wall_plane p75(array[count][1],array[count][2],array[count][3],array[count][0]);
				con.add_wall(p75);
				continue;
			}
			if(count==76)
			{
				wall_plane p76(array[count][1],array[count][2],array[count][3],array[count][0]);
				con.add_wall(p76);
				continue;
			}
			if(count==77)
			{
				wall_plane p77(array[count][1],array[count][2],array[count][3],array[count][0]);
				con.add_wall(p77);
				continue;
			}
			if(count==78)
			{
				wall_plane p78(array[count][1],array[count][2],array[count][3],array[count][0]);
				con.add_wall(p78);
				continue;
			}
			if(count==79)
			{
				wall_plane p79(array[count][1],array[count][2],array[count][3],array[count][0]);
				con.add_wall(p79);
				continue;
			}
			if(count==80)
			{
				wall_plane p80(array[count][1],array[count][2],array[count][3],array[count][0]);
				con.add_wall(p80);
				continue;
			}
			if(count==81)
			{
				wall_plane p81(array[count][1],array[count][2],array[count][3],array[count][0]);
				con.add_wall(p81);
				continue;
			}
			if(count==82)
			{
				wall_plane p82(array[count][1],array[count][2],array[count][3],array[count][0]);
				con.add_wall(p82);
				continue;
			}
			if(count==83)
			{
				wall_plane p83(array[count][1],array[count][2],array[count][3],array[count][0]);
				con.add_wall(p83);
				continue;
			}
			if(count==84)
			{
				wall_plane p84(array[count][1],array[count][2],array[count][3],array[count][0]);
				con.add_wall(p84);
				continue;
			}
			if(count==85)
			{
				wall_plane p85(array[count][1],array[count][2],array[count][3],array[count][0]);
				con.add_wall(p85);
				continue;
			}
			if(count==86)
			{
				wall_plane p86(array[count][1],array[count][2],array[count][3],array[count][0]);
				con.add_wall(p86);
				continue;
			}
			if(count==87)
			{
				wall_plane p87(array[count][1],array[count][2],array[count][3],array[count][0]);
				con.add_wall(p87);
				continue;
			}
			if(count==88)
			{
				wall_plane p88(array[count][1],array[count][2],array[count][3],array[count][0]);
				con.add_wall(p88);
				continue;
			}
			if(count==89)
			{
				wall_plane p89(array[count][1],array[count][2],array[count][3],array[count][0]);
				con.add_wall(p89);
				continue;
			}
			if(count==90)
			{
				wall_plane p90(array[count][1],array[count][2],array[count][3],array[count][0]);
				con.add_wall(p90);
				continue;
			}
			if(count==91)
			{
				wall_plane p91(array[count][1],array[count][2],array[count][3],array[count][0]);
				con.add_wall(p91);
				continue;
			}
			if(count==92)
			{
				wall_plane p92(array[count][1],array[count][2],array[count][3],array[count][0]);
				con.add_wall(p92);
				continue;
			}
			if(count==93)
			{
				wall_plane p93(array[count][1],array[count][2],array[count][3],array[count][0]);
				con.add_wall(p93);
				continue;
			}
			if(count==94)
			{
				wall_plane p94(array[count][1],array[count][2],array[count][3],array[count][0]);
				con.add_wall(p94);
				continue;
			}
			if(count==95)
			{
				wall_plane p95(array[count][1],array[count][2],array[count][3],array[count][0]);
				con.add_wall(p95);
				continue;
			}
			if(count==96)
			{
				wall_plane p96(array[count][1],array[count][2],array[count][3],array[count][0]);
				con.add_wall(p96);
				continue;
			}
			if(count==97)
			{
				wall_plane p97(array[count][1],array[count][2],array[count][3],array[count][0]);
				con.add_wall(p97);
				continue;
			}
			if(count==98)
			{
				wall_plane p98(array[count][1],array[count][2],array[count][3],array[count][0]);
				con.add_wall(p98);
				continue;
			}
			if(count==99)
			{
				wall_plane p99(array[count][1],array[count][2],array[count][3],array[count][0]);
				con.add_wall(p99);
				continue;
			}
			if(count==100)
			{
				wall_plane p100(array[count][1],array[count][2],array[count][3],array[count][0]);
				con.add_wall(p100);
				continue;
			}
			if(count==101)
			{
				wall_plane p101(array[count][1],array[count][2],array[count][3],array[count][0]);
				con.add_wall(p101);
				continue;
			}
			if(count==102)
			{
				wall_plane p102(array[count][1],array[count][2],array[count][3],array[count][0]);
				con.add_wall(p102);
				continue;
			}
			if(count==103)
			{
				wall_plane p103(array[count][1],array[count][2],array[count][3],array[count][0]);
				con.add_wall(p103);
				continue;
			}
			if(count==104)
			{
				wall_plane p104(array[count][1],array[count][2],array[count][3],array[count][0]);
				con.add_wall(p104);
				continue;
			}
			if(count==105)
			{
				wall_plane p105(array[count][1],array[count][2],array[count][3],array[count][0]);
				con.add_wall(p105);
				continue;
			}
			if(count==106)
			{
				wall_plane p106(array[count][1],array[count][2],array[count][3],array[count][0]);
				con.add_wall(p106);
				continue;
			}
			if(count==107)
			{
				wall_plane p107(array[count][1],array[count][2],array[count][3],array[count][0]);
				con.add_wall(p107);
				continue;
			}
			if(count==108)
			{
				wall_plane p108(array[count][1],array[count][2],array[count][3],array[count][0]);
				con.add_wall(p108);
				continue;
			}
			if(count==109)
			{
				wall_plane p109(array[count][1],array[count][2],array[count][3],array[count][0]);
				con.add_wall(p109);
				continue;
			}
			if(count==110)
			{
				wall_plane p110(array[count][1],array[count][2],array[count][3],array[count][0]);
				con.add_wall(p110);
				continue;
			}
			if(count==111)
			{
				wall_plane p111(array[count][1],array[count][2],array[count][3],array[count][0]);
				con.add_wall(p111);
				continue;
			}
			if(count==112)
			{
				wall_plane p112(array[count][1],array[count][2],array[count][3],array[count][0]);
				con.add_wall(p112);
				continue;
			}
			if(count==113)
			{
				wall_plane p113(array[count][1],array[count][2],array[count][3],array[count][0]);
				con.add_wall(p113);
				continue;
			}
			if(count==114)
			{
				wall_plane p114(array[count][1],array[count][2],array[count][3],array[count][0]);
				con.add_wall(p114);
				continue;
			}
			if(count==115)
			{
				wall_plane p115(array[count][1],array[count][2],array[count][3],array[count][0]);
				con.add_wall(p115);
				continue;
			}
			if(count==116)
			{
				wall_plane p116(array[count][1],array[count][2],array[count][3],array[count][0]);
				con.add_wall(p116);
				continue;
			}
			if(count==117)
			{
				wall_plane p117(array[count][1],array[count][2],array[count][3],array[count][0]);
				con.add_wall(p117);
				continue;
			}
			if(count==118)
			{
				wall_plane p118(array[count][1],array[count][2],array[count][3],array[count][0]);
				con.add_wall(p118);
				continue;
			}
			if(count==119)
			{
				wall_plane p119(array[count][1],array[count][2],array[count][3],array[count][0]);
				con.add_wall(p119);
				continue;
			}
			if(count==120)
			{
				wall_plane p120(array[count][1],array[count][2],array[count][3],array[count][0]);
				con.add_wall(p120);
				continue;
			}
			if(count==121)
			{
				wall_plane p121(array[count][1],array[count][2],array[count][3],array[count][0]);
				con.add_wall(p121);
				continue;
			}
			if(count==122)
			{
				wall_plane p122(array[count][1],array[count][2],array[count][3],array[count][0]);
				con.add_wall(p122);
				continue;
			}
			if(count==123)
			{
				wall_plane p123(array[count][1],array[count][2],array[count][3],array[count][0]);
				con.add_wall(p123);
				continue;
			}
			if(count==124)
			{
				wall_plane p124(array[count][1],array[count][2],array[count][3],array[count][0]);
				con.add_wall(p124);
				continue;
			}
			if(count==125)
			{
				wall_plane p125(array[count][1],array[count][2],array[count][3],array[count][0]);
				con.add_wall(p125);
				continue;
			}
			if(count==126)
			{
				wall_plane p126(array[count][1],array[count][2],array[count][3],array[count][0]);
				con.add_wall(p126);
				continue;
			}
			if(count==127)
			{
				wall_plane p127(array[count][1],array[count][2],array[count][3],array[count][0]);
				con.add_wall(p127);
				continue;
			}
			if(count==128)
			{
				wall_plane p128(array[count][1],array[count][2],array[count][3],array[count][0]);
				con.add_wall(p128);
				continue;
			}
			if(count==129)
			{
				wall_plane p129(array[count][1],array[count][2],array[count][3],array[count][0]);
				con.add_wall(p129);
				continue;
			}
			if(count==130)
			{
				wall_plane p130(array[count][1],array[count][2],array[count][3],array[count][0]);
				con.add_wall(p130);
				continue;
			}
		}
		
		// Variables for the tessellation
		int i=0;
		double x,y,z;
		
		// Inserting the Voronoi points in the polyhedron in a conditioned method
		if (a_voro<=1){
			double l_box;
			int n_elements;
			
			// Insert particles into the container in a tetahedric grid, checking  if ...
			// efectively the point lies inside the polyhedron
			
			l_box = max(max((x_max - x_min), (y_max - y_min)),(z_max - z_min));
			n_elements = ceil(pow(n_sub_par,(0.333333333)));
			
			l_box = l_box/n_elements;
			
			// Coordinates of the voronoi points
			for (int count_z=0; count_z<=n_elements; count_z++){
				z = z_min + l_box/2 + l_box*count_z + l_box*a_voro*(rnd()-0.5);
				for (int count_y=0; count_y<=n_elements; count_y++){
					y = y_min + l_box/2 + l_box*count_y + l_box*a_voro*(rnd()-0.5);
					for (int count_x=0; count_x<=n_elements; count_x++){
						x = x_min + l_box/2 + l_box*count_x + l_box*a_voro*(rnd()-0.5);
						// Checking if the point is inside the body
						if (con.point_inside(x,y,z)){
							con.put(i,x,y,z);
							i++;
						}
					}
				}
			}
		}
		else if(a_voro==2){
			while(i<n_sub_par) {
				x=x_min+rnd()*(x_max-x_min);
				y=y_min+rnd()*(y_max-y_min);
				z=z_min+rnd()*(z_max-z_min);
				if (con.point_inside(x,y,z)) {
					con.put(i,x,y,z);
					i++;
				}
			}
		}
		else {
			cout << "Unknown tessellation procedure" << "\n";
			return 0;
		}
		
		// Writing in the terminal
		cout << "########################################" << "\n";
		cout << "         Succesful tessellation         " << "\n";
		
		// Output the particle positions and the Voronoi cells in POV-Ray format
		
		// Writing file for POVRAY that containes the voronoi points information
		string name_file_p = "./POSTPRO/VORONOI/pray_points" + to_string(parsing) +".pov";
		FILE *fp = safe_fopen(&name_file_p[0], "w");
		con.draw_particles_pov(fp);
		fclose(fp);
		
		// Writing file for POVRAY that contains the vertices and edges information
		string name_file_v = "./POSTPRO/VORONOI/pray_vertices" + to_string(parsing) +".pov";
		FILE *fv = safe_fopen(&name_file_v[0], "w");
		con.draw_cells_pov(fv);
		fclose(fv);
		
		// Custom file to be read by the python script
		// Contains the sub-particle ID, centroid, # of vertices, and vertices cordinates
		string name_file_custom = "./POSTPRO/VORONOI/vertices" + to_string(parsing);
		FILE *fc = safe_fopen(&name_file_custom[0], "w");
		con.print_custom("ID=%i, Center=%C, N_vertices=%w, Vertices=%P", fc);
		fclose(fc);
		
		// Custom file for the post processing of the tessellation
		// Contains the particle ID, # of vertices, # of faces, # of edges, and the volume.
		string name_file_postpro1 = "./POSTPRO/VORONOI/post_pro" + to_string(parsing);
		FILE *fpo1 = safe_fopen(&name_file_postpro1[0], "w");
		con.print_custom("%i, %w, %s, %g, %v", fpo1);
		fclose(fpo1);
		
		// Custom file for the post processing of the tessellation
		// Contains the neighbors of each sub-particle
		string name_file_postpro2 = "./POSTPRO/VORONOI/post_pro_nei" + to_string(parsing);
		FILE *fpo2 = safe_fopen(&name_file_postpro2[0], "w");
		con.print_custom("%n", fpo2);
		fclose(fpo2);
		
		// Custom file for the post processing of the tessellation
		// Contains the order of the faces, and the area of the faces
		string name_file_postpro3 = "./POSTPRO/VORONOI/post_pro_afaces" + to_string(parsing);
		FILE *fpo3 = safe_fopen(&name_file_postpro3[0], "w");
		con.print_custom("%a, %f", fpo3);
		fclose(fpo3);
		
		
	}
}
