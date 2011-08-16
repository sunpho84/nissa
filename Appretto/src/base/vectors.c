#pragma once

//check weather the borders are allocated
void check_minimal_allocated_size(void *v,int l,const char *err_mess)
{
  appretto_vect *vect=v-sizeof(appretto_vect);
  if(vect->nel<l)
    if(rank==0)
      {
	fprintf(stderr,"Error \"%s\" in ",err_mess);
	appretto_vect_content_fprintf(stderr,vect);
      }
}

//check weather the borders are allocated
void check_borders_allocated(void *v)
{check_minimal_allocated_size(v,loc_vol+loc_bord,"border not allocated");}

//check weather the borders and edges are allocated
void check_edges_allocated(void *v)
{check_minimal_allocated_size(v,loc_vol+loc_bord+loc_edge,"border not allocated");}

//print the content of an appretto vect
void appretto_vect_content_fprintf(FILE *fout,appretto_vect *vect)
{
  if(rank==0)
    {
      fprintf(fout,"\"%s\" ",vect->tag);
      fprintf(fout,"of %d elements of type \"%s\" (%d bytes) ",vect->nel,vect->type,vect->nel*vect->size_per_el);
      fprintf(fout,"allocated in file %s line %d\n",vect->file,vect->line);
    }
}
//wrappers
void appretto_vect_content_printf(appretto_vect *vect)
{appretto_vect_content_fprintf(stdout,vect);}
void last_appretto_vect_content_printf()
{
  if(rank==0) printf("Last appretto vect content: ");
  appretto_vect_content_printf(last_appretto_vect);
}

//print all appretto vect
void print_all_appretto_vect_content()
{
  appretto_vect *curr=&(main_appretto_vect);
  do
    {  
      appretto_vect_content_printf(curr);
      curr=curr->next;
    }
  while(curr!=NULL);
}

//print all appretto vect
int compute_appretto_vect_memory_usage()
{
  int tot=0;
  appretto_vect *curr=&(main_appretto_vect);
  do
    {  
      tot+=curr->nel*curr->size_per_el;
      curr=curr->next;
    }
  while(curr!=NULL);
  
  return tot;
}

//initialize the first vector
void initialize_main_appretto_vect()
{
  if_appretto_vect_not_initialized()
    {
      appretto_required_memory=0;
      appretto_max_required_memory=0;
      last_appretto_vect=&main_appretto_vect;
      sprintf(main_appretto_vect.tag,"base");
      sprintf(main_appretto_vect.type,"(null)");
      main_appretto_vect.prev=main_appretto_vect.next=NULL;
      main_appretto_vect.nel=0;
      main_appretto_vect.size_per_el=0;
      memcpy(main_appretto_vect.file,__FILE__+max_int(0,strlen(__FILE__)-12),12);
      main_appretto_vect.line=__LINE__;
      main_appretto_arr=(void*)last_appretto_vect+sizeof(appretto_vect);
    }
}

//allocate an appretto vector 
void *appretto_true_malloc(const char *tag,int nel,int size_per_el,const char *type,const char *file,int line)
{
  if_appretto_vect_not_initialized() initialize_main_appretto_vect();
  
  int size=nel*size_per_el;
  //try to allocate the new vector
  appretto_vect *new=(void*)malloc(size+sizeof(appretto_vect));
  if(new==NULL)
    if(rank==0)
      {
        fprintf(stderr,"Error! Could not allocate vector named \"%s\" of %d elements of types %s\n",tag,nel,type);
        fprintf(stderr,"(total size: %d bytes) request on line %d of file %s\n",size,line,file);
	MPI_Abort(MPI_COMM_WORLD,1);
      }
  
  //fill the vector with information supplied
  new->line=line;
  new->nel=nel;
  new->size_per_el=size_per_el;
  take_last_characters(new->file,file,appretto_vect_string_length);
  take_last_characters(new->tag,tag,appretto_vect_string_length);
  take_last_characters(new->type,type,appretto_vect_string_length);
  
  //append the new vector to the list
  new->next=NULL;
  new->prev=last_appretto_vect;
  
  last_appretto_vect->next=new;
  last_appretto_vect=new;
  
  if(rank==0 && debug>1) 
    {
      printf("Allocated vector ");
      appretto_vect_content_printf(last_appretto_vect);
    }
  
  //Update the amount of required memory
  appretto_required_memory+=size;
  appretto_max_required_memory=max_int(appretto_max_required_memory,appretto_required_memory);
  
  return (void*)last_appretto_vect+sizeof(appretto_vect);
}

//release a vector
void* appretto_true_free(void *arr,const char *file,int line)
{
  if(arr!=NULL)
    {
      appretto_vect *vect=(appretto_vect*)(arr-sizeof(appretto_vect));
      appretto_vect *prev=vect->prev;
      appretto_vect *next=vect->next;
      
      if(rank==0 && debug>1)
	{
	  printf("At line %d of file %s freeing vector ",line,file);
	  appretto_vect_content_printf(vect);
	}
  
      //detach from previous
      prev->next=next;
      
      //if not last element
      if(next!=NULL) next->prev=prev;
      else last_appretto_vect=prev;
    
      //update the appretto required memory
      appretto_required_memory-=(vect->size_per_el*vect->nel);
  
      free(vect);
    }
  else crash("Error, trying to delocate a NULL vector on line: %d of file: %s\n",line,file);
  
  return NULL;
}
