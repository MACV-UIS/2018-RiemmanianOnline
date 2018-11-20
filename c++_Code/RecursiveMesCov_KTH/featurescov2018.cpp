#include <iostream>
#include <vector>
#include <sstream>
#include <fstream>
#include <math.h>
#include <stdlib.h>

using namespace std;

string ToString(int num)
{
    stringstream stream;
    stream << num;
    return stream.str();
}

string ToStringf(float num)
{
    stringstream stream;
    stream << num;
    return stream.str();
}

vector<string> splitStr(string str, char delimiter)
{
    vector<string> internal;
    stringstream ss(str);
    string tok;
    while(getline(ss, tok, delimiter)){
        internal.push_back(tok);
    }
    return internal;
}


int main(int argc, char** argv){

//***********************************************
//float*** coordenadas;
//int x,y,z,m;
int cont_tj, cont_pt_tj, cont_z;
int numTotalFrames, heightFrame, widthFrame, num_tot_traj;
//int inicio_traj, fin_traj;
//int frame_num;
//***********************************************

vector<string> vector_metadatos,vector_auxiliar;
vector<string> vector_metadatos2;

vector<string> vector_ftconcat;
vector<string> vector_ftconcat_norm;
//vector<string> vector_ftconcat_wonly;
//vector<string> vector_ftconcat_normonly;
vector<double> vector_ftconcat_num;

vector<string> vector_valores, vector_valores2;
vector< vector<string> > vector_total, vector_total2;

vector<string> vector_valores_selected;
vector< vector<string> > vector_total_selected;

vector<string> vectortemp;
vector< vector<string> > vectorcov, vectorcov2;

string read_line;

string ruta = argv[1];      //current test file
string myfilename;
string myfilenamenoext;

string ruta2;
string ruta3;

string acumft,acumft2;

vector_metadatos2 = splitStr(ruta,'/');
myfilename = vector_metadatos2[9];
vector_metadatos2.resize(0);
vector_metadatos2.clear();

vector_metadatos2 = splitStr(myfilename,'.');
myfilenamenoext = vector_metadatos2[0];
vector_metadatos2.resize(0);
vector_metadatos2.clear();

cout<<"ruta: "<<ruta<<endl;
cout<<"myfilename: "<<myfilename<<endl;
cout<<"myfilenamenoext: "<<myfilenamenoext<<endl;


//string rutatraj = argv[2];

//string rutatrajfile = argv[2]+myfilename;
//cout<<"rutatrajfile: "<<rutatrajfile<<endl;
//ruta3 = argv[3]+myfilename;


int t;
int limite;

ifstream datafiles_stream2(ruta.c_str(), fstream::in);
while(datafiles_stream2.good()){
	getline(datafiles_stream2, read_line);
	if(read_line.size()>0){

		vector_metadatos = splitStr(read_line,' ');
		limite = (vector_metadatos.size()-1)/2;

		for(t=1; t<=limite; t++){
			vector_metadatos2 = splitStr(vector_metadatos[t],':');
			vector_valores.push_back(vector_metadatos2[1]);
			vector_metadatos2.clear();
		}
		vector_total.push_back(vector_valores);
		vector_valores.clear();

		for(t=limite+1; t<vector_metadatos.size(); t++){
			vector_metadatos2 = splitStr(vector_metadatos[t],':');
			vector_valores2.push_back(vector_metadatos2[1]);
			vector_metadatos2.clear();
		}
		vector_total2.push_back(vector_valores2);
		vector_valores2.clear();

		vector_metadatos.clear();

	}
}
datafiles_stream2.close();




int NFEATURES = 12;


for(int i=0; i< NFEATURES; i++){
	for(int j=0; j< NFEATURES; j++){
		vectortemp.push_back("");
	}
	vectorcov.push_back(vectortemp);
}
vectortemp.clear();

for(int i=0; i< NFEATURES; i++){
	for(int j=0; j< NFEATURES; j++){
		vectortemp.push_back("");
	}
	vectorcov2.push_back(vectortemp);
}
vectortemp.clear();







int celda = 0;
int selected = 0;
int zeros = 0;
int fila = 0;

for(t=0; t<vector_total.size(); t++){
	for(celda=0; celda<vector_total[0].size(); t++){
		if(zeros!=0){
			for(int z=0; z<zeros; z++){
				vectorcov[fila][z] = "0";
			}
		}
		vector_total[t][celda]

for(int i=0; i< NFEATURES; i++){
	if(i != selected){
		for(int j=i; j< NFEATURES; j++){
			if(j != selected){
				vector_valores_selected.push_back(vector_total[t][celda]);
			}
		}
	}
	celda++;
}

vector_total_selected.push_back(vector_valores_selected);
vector_valores_selected.clear();

	}
}
















int cont_P;

string ruta4 = "new-"+myfilenamenoext+".txt";
ofstream outdata_trajF(ruta4.c_str(), fstream::out);

for(t=0; t<vector_total_selected.size(); t++){

	cont_P = 1;
	for(int i=0; i< vector_total_selected[0].size(); i++){
		outdata_trajF<<cont_P<<":"<<vector_total_selected[t][i]<<" ";
		cont_P++;
	}
	outdata_trajF<<endl;

}
outdata_trajF.close();






return 0;

}
