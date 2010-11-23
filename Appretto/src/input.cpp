#pragma once

ifstream input_global;

void open_input(char *input_path)
{
  if(rank==0)
    {
      input_global.open(input_path);
      if(input_global.good()!=1)
	{
	  cerr<<"File '"<<input_path<<"' not found"<<endl;
	  MPI_Abort(MPI_COMM_WORLD,1);
	}

      if(debug) cout<<"File '"<<input_path<<"' opened"<<endl;
    }	
}

void close_input()
{
  if(rank==0) input_global.close();
}

//Read an integer from the file
void read_int(int &in)
{
  if(rank==0) input_global>>in;
  MPI_Bcast(&in,1,MPI_INT,0,MPI_COMM_WORLD);
}

//Read a string from the file
void read_str(char *str,int length)
{
  if(rank==0) input_global>>str;
  MPI_Bcast(str,length,MPI_BYTE,0,MPI_COMM_WORLD);
}

//Read a string from the file and check against the argument
void expect_str(const char *exp_str)
{
  char obt_str[1024];

  read_str(obt_str,1024);
  
  if(strcasecmp(exp_str,obt_str)!=0)
    {
      if(rank==0) cerr<<"Error, expexcted '"<<exp_str<<"' in input file, obtained: '"<<obt_str<<"'"<<endl;
      MPI_Abort(MPI_COMM_WORLD,1);
    }
}

//Read an integer checking the tag
void read_int(const char *exp_str,int &in)
{
  expect_str(exp_str);
  read_int(in);

  if(rank==0) cout<<"Read variable '"<<exp_str<<"' with value: "<<in<<endl;
}

//Read a string checking the tag
void read_str(const char *exp_str,char *in,int length=1024)
{
  expect_str(exp_str);
  read_str(in,length);

  if(rank==0) cout<<"Read variable '"<<exp_str<<"' with value: "<<in<<endl;
}
