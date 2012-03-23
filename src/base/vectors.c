#pragma once

//check weather the borders are allocated
void check_minimal_allocated_size(void *v,int l,const char *err_mess)
{
  nissa_vect *vect=(nissa_vect*)((char*)v-sizeof(nissa_vect));
  if(vect->nel<l)
    if(rank==0)
      {
	fprintf(stderr,"Error \"%s\" in ",err_mess);
	nissa_vect_content_fprintf(stderr,vect);
      }
}

//check weather the borders are allocated
void check_lx_borders_allocated(void *v)
{check_minimal_allocated_size(v,loc_vol+loc_bord,"lx border not allocated");}
void check_eo_borders_allocated(void *v)
{check_minimal_allocated_size(v,loc_volh+loc_bordh,"eo border not allocated");}

//check weather the borders and edges are allocated
void check_lx_edges_allocated(void *v)
{check_minimal_allocated_size(v,loc_vol+loc_bord+loc_edge,"lx edges not allocated");}

//return the name of the vector
char *get_vect_name(void *v)
{
  nissa_vect *vect=(nissa_vect*)((char*)v-sizeof(nissa_vect));
  return vect->tag;
}

//print the content of an nissa vect
void nissa_vect_content_fprintf(FILE *fout,nissa_vect *vect)
{
  if(rank==0)
    {
      fprintf(fout,"\"%s\" ",vect->tag);
      fprintf(fout,"of %d elements of type \"%s\" (%d bytes) ",vect->nel,vect->type,vect->nel*vect->size_per_el);
      fprintf(fout,"allocated in file %s line %d\n",vect->file,vect->line);
    }
}
//wrappers
void nissa_vect_content_printf(nissa_vect *vect)
{nissa_vect_content_fprintf(stdout,vect);}
void last_nissa_vect_content_printf()
{
  master_printf("Last nissa vect content: ");
  nissa_vect_content_printf(last_nissa_vect);
}

//print all nissa vect
void print_all_nissa_vect_content()
{
  nissa_vect *curr=&(main_nissa_vect);
  do
    {  
      nissa_vect_content_printf(curr);
      curr=curr->next;
    }
  while(curr!=NULL);
}

//print all nissa vect
int compute_nissa_vect_memory_usage()
{
  int tot=0;
  nissa_vect *curr=&(main_nissa_vect);
  do
    {  
      tot+=curr->nel*curr->size_per_el;
      curr=curr->next;
    }
  while(curr!=NULL);
  
  return tot;
}

//initialize the first vector
void initialize_main_nissa_vect()
{
  if_nissa_vect_not_initialized()
    {
      nissa_required_memory=0;
      nissa_max_required_memory=0;
      last_nissa_vect=&main_nissa_vect;
      sprintf(main_nissa_vect.tag,"base");
      sprintf(main_nissa_vect.type,"(null)");
      main_nissa_vect.prev=main_nissa_vect.next=NULL;
      main_nissa_vect.nel=0;
      main_nissa_vect.size_per_el=0;
      memcpy(main_nissa_vect.file,__FILE__+max_int(0,strlen(__FILE__)-12),12);
      main_nissa_vect.line=__LINE__;
      main_nissa_arr=(char*)last_nissa_vect+sizeof(nissa_vect);

      master_printf("Vector memory manager started.\n");
    }
}

//allocate an nissa vector 
void *nissa_true_malloc(const char *tag,int nel,int size_per_el,const char *type,const char *file,int line)
{
  if_nissa_vect_not_initialized() initialize_main_nissa_vect();
  
  int size=nel*size_per_el;
  //try to allocate the new vector
  nissa_vect *new=(void*)malloc(size+sizeof(nissa_vect));
  if(new==NULL)
    crash("could not allocate vector named \"%s\" of %d elements of type %s (total size: %d bytes) request on line %d of file %s"
	  ,                                  tag,    nel,                type,           size,                     line,      file);
  
  //fill the vector with information supplied
  new->line=line;
  new->nel=nel;
  new->size_per_el=size_per_el;
  new->flag=0;
  take_last_characters(new->file,file,nissa_vect_string_length);
  take_last_characters(new->tag,tag,nissa_vect_string_length);
  take_last_characters(new->type,type,nissa_vect_string_length);
  
  //append the new vector to the list
  new->next=NULL;
  new->prev=last_nissa_vect;
  
  last_nissa_vect->next=new;
  last_nissa_vect=new;
  
  if(debug_lvl>1) 
    {
      master_printf("Allocated vector ");
      nissa_vect_content_printf(last_nissa_vect);
    }
  
  //define returned pointer and check for its alignement
  void *return_ptr=(void*)last_nissa_vect+sizeof(nissa_vect);
  int offset=((long long int)(return_ptr))%nissa_vect_alignment;
  if(offset!=0)
    crash("memory alignment problem, vector %s has %d offset",tag,offset);
  
  //Update the amount of required memory
  nissa_required_memory+=size;
  nissa_max_required_memory=max_int(nissa_max_required_memory,nissa_required_memory);
  
  return return_ptr;
}

//release a vector
void* nissa_true_free(void *arr,const char *file,int line)
{
  if(arr!=NULL)
    {
      nissa_vect *vect=(nissa_vect*)((char*)arr-sizeof(nissa_vect));
      nissa_vect *prev=vect->prev;
      nissa_vect *next=vect->next;
      
      if(debug_lvl>1)
	{
	  master_printf("At line %d of file %s freeing vector ",line,file);
	  nissa_vect_content_printf(vect);
	}
  
      //detach from previous
      prev->next=next;
      
      //if not last element
      if(next!=NULL) next->prev=prev;
      else last_nissa_vect=prev;
    
      //update the nissa required memory
      nissa_required_memory-=(vect->size_per_el*vect->nel);
  
      free(vect);
    }
  else crash("Error, trying to delocate a NULL vector on line: %d of file: %s\n",line,file);
  
  return NULL;
}
