#include <include.h>

int T;
int njack=16;

int *icombo;
char **data_path;
int tmin,tmax;

int ndata;

jvec load_2pts(int idata)
{
  return jvec_load(data_path[idata],T,njack,icombo[idata]);
}

void read_input()
{
  FILE *input_file=open_file("input","r");

  read_formatted_from_file_expecting((char*)(&ndata),input_file,"%d","ndata");
  data_path=(char**)malloc(ndata*sizeof(char*));
  icombo=(int*)malloc(ndata*sizeof(int));
  for(int idata=0;idata<ndata;idata++)
    {
      data_path[idata]=(char*)malloc(1024);
      read_formatted_from_file_expecting(data_path[idata],input_file,"%s","data_path");
      read_formatted_from_file((char*)&(icombo[idata]),input_file,"%d","icombo");
    }
  
  read_formatted_from_file_expecting((char*)(&T),input_file,"%d","T");
  
  read_formatted_from_file_expecting((char*)(&tmin),input_file,"%d","tmin");
  read_formatted_from_file_expecting((char*)(&tmax),input_file,"%d","tmax");
  
  fclose(input_file);
}

int main()
{
  read_input();
  
  ///////////////////////////// Load two points for standing D and D* //////////////////////////
  
  //load sl
  jvec C(T,njack);
  
  for(int idata=0;idata<ndata;idata++)
    {
      C=(C*idata+load_2pts(idata))/(idata+1);
      
      //compute mass
      jack M=constant_fit(effective_mass(C.simmetrized(1)),tmin,tmax,combine("M_%02d.xmg",idata).c_str());
      cout<<"nsources: "<<idata<<" mass: "<<M<<endl;
      
      if(idata==ndata-1) M.write_to_binfile("M");
    }
  
  return 0;
}
