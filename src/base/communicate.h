#ifndef _COMMUNICATE_H
#define _COMMUNICATE_H

#include <mpi.h>

#include "../new_types/new_types_definitions.h"

void communicate_eo_borders(char **data,MPI_Datatype *MPI_EO_BORDS_SEND_TXY,MPI_Datatype *MPI_EV_BORDS_SEND_Z,MPI_Datatype *MPI_OD_BORDS_SEND_Z,MPI_Datatype *MPI_EO_BORDS_RECE,int nbytes_per_site);
void communicate_eo_color_borders(color **eos);
void communicate_eo_edges(char **data,MPI_Datatype *MPI_EO_BORDS_SEND_TXY,MPI_Datatype *MPI_EV_BORDS_SEND_Z,MPI_Datatype *MPI_OD_BORDS_SEND_Z,MPI_Datatype *MPI_EO_BORDS_RECE,MPI_Datatype *MPI_EDGES_SEND,MPI_Datatype *MPI_EDGES_RECE,int nbytes_per_site);
void communicate_eo_quad_su3_borders(quad_su3 **eo_conf);
void communicate_eo_quad_su3_edges(quad_su3 **conf);
void communicate_ev_borders(char *ev_data,MPI_Datatype *MPI_EV_BORDS_SEND_TXY,MPI_Datatype *MPI_EV_BORDS_SEND_Z,MPI_Datatype *MPI_EV_BORDS_RECE,int nbytes_per_site);
void communicate_ev_color_borders(color *ev);
void communicate_ev_spin_borders(spin *ev);
void communicate_ev_spincolor_borders(spincolor *ev);
void communicate_lx_borders(char *data,MPI_Datatype *MPI_BORDS_SEND,MPI_Datatype *MPI_BORDS_RECE,int nbytes_per_site);
void communicate_lx_edges(char *data,MPI_Datatype *MPI_BORDS_SEND,MPI_Datatype *MPI_BORDS_RECE,MPI_Datatype *MPI_EDGES_SEND,MPI_Datatype *MPI_EDGES_RECE,int nbytes_per_site);
void communicate_lx_quad_su3_borders(quad_su3 *conf);
void communicate_lx_quad_su3_edges(quad_su3 *conf);
void communicate_lx_spincolor_128_borders(spincolor_128 *s);
void communicate_lx_spin_borders(spin *s);
void communicate_lx_color_borders(color *s);
void communicate_lx_spinspin_borders(spinspin *s);
void communicate_lx_spincolor_borders(spincolor *s);
void communicate_lx_su3_borders(su3 *u);
void communicate_lx_su3_edges(su3 *u);
void communicate_od_borders(char *od_data,MPI_Datatype *MPI_EV_BORDS_SEND_TXY,MPI_Datatype *MPI_OD_BORDS_SEND_Z,MPI_Datatype *MPI_EV_BORDS_RECE,int nbytes_per_site);
void communicate_od_color_borders(color *od);
void communicate_od_spin_borders(spin *od);
void communicate_od_spincolor_borders(spincolor *od);
void communicate_ev_spincolor_128_borders(spincolor_128 *ev);
void communicate_od_spincolor_128_borders(spincolor_128 *od);
void finish_communicating_ev_borders(int &nrequest,MPI_Request *request,char *ev_data);
void finish_communicating_ev_color_borders(int &nrequest,MPI_Request *request,color *ev);
void finish_communicating_ev_spincolor_borders(int &nrequest,MPI_Request *request,spincolor *ev);
void finish_communicating_lx_borders(int &nrequest,MPI_Request *request,char *data);
void finish_communicating_lx_spincolor_borders(int &nrequest,MPI_Request *request,spincolor *s);
void start_communicating_ev_borders(int &nrequest_ret,MPI_Request *request,char *ev_data,MPI_Datatype *MPI_EV_BORDS_SEND_TXY,MPI_Datatype *MPI_EV_BORDS_SEND_Z,MPI_Datatype *MPI_EV_BORDS_RECE,int nbytes_per_site);
void start_communicating_ev_color_borders(int &nrequest,MPI_Request *request,color *ev);
void start_communicating_ev_spincolor_borders(int &nrequest,MPI_Request *request,spincolor *ev);
void start_communicating_lx_borders(int &nrequest_ret,MPI_Request *request,char *data,MPI_Datatype *MPI_BORDS_SEND,MPI_Datatype *MPI_BORDS_RECE,int nbytes_per_site);
void start_communicating_lx_spincolor_borders(int &nrequest,MPI_Request *request,spincolor *s);

#endif
