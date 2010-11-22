#pragma once

ifstream input_global;

void open_input(const char *input_path)
{
  if(rank==0) input_global.open(input_path);
}

void close_input()
{
  if(rank==0) input_global.close();
}

void read_int(int &in)
{
  if(rank==0) input_global>>in;
  MPI_Bcast(&in,1,MPI_INT,0,MPI_COMM_WORLD);
}

void read_int(const char *exp_str,int &in)
{
  char obt_str[1024];

  if(rank==0) 
    {
      input_global>>obt_str;

      if(strcmp(exp_str,obt_str)!=0)
	{
	  cerr<<"Error, expexcted '"<<exp_str<<"' in input file, obtained: '"<<obt_str<<"'"<<endl;
	  MPI_Abort(MPI_COMM_WORLD,1);
	}
    }
  read_int(in);

  if(rank==0) cout<<"Read variable '"<<exp_str<<"' with value: "<<in<<endl;
}
