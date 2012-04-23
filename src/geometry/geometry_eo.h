#pragma once

void addrem_stagphases_to_eo_conf(quad_su3 **eo_conf);
void initialize_eo_bord_receivers_of_kind(MPI_Datatype *MPI_EO_BORD_RECE,MPI_Datatype *base);
void initialize_eo_bord_senders_of_kind(MPI_Datatype *MPI_EO_BORD_SEND_TXY,MPI_Datatype *MPI_EV_BORD_SEND_Z,MPI_Datatype *MPI_OD_BORD_SEND_Z,MPI_Datatype *base);
void initialize_eo_edge_receivers_of_kind(MPI_Datatype *MPI_EDGES_RECE,MPI_Datatype *base);
void initialize_eo_edge_senders_of_kind(MPI_Datatype *MPI_EO_EDGES_SEND,MPI_Datatype *base);
void paste_eo_parts_into_lx_color(color *out_lx,color **in_eo);
void paste_eo_parts_into_lx_conf(quad_su3 *out_lx,quad_su3 **in_eo);
void paste_eo_parts_into_lx_spincolor(spincolor *out_lx,spincolor **in_eo);
void paste_eo_parts_into_lx_vector(char *out_lx,char **in_eo,int bps);
void set_eo_bord_senders_and_receivers(MPI_Datatype *MPI_EO_BORD_SEND_TXY,MPI_Datatype *MPI_EV_BORD_SEND_Z,MPI_Datatype *MPI_OD_BORD_SEND_Z,MPI_Datatype *MPI_EO_BORD_RECE,MPI_Datatype *base);
void set_eo_edge_senders_and_receivers(MPI_Datatype *MPI_EO_EDGES_SEND,MPI_Datatype *MPI_EO_EDGES_RECE,MPI_Datatype *base);
void set_eo_geometry();
void split_lx_color_into_eo_parts(color **eo_out,color *lx_in);
void split_lx_conf_into_eo_parts(quad_su3 **eo_out,quad_su3 *lx_in);
void split_lx_spincolor_into_eo_parts(spincolor **eo_out,spincolor *lx_in);
void split_lx_vector_into_eo_parts(char **out_eo,char *in_lx,int bps);
void take_e_or_o_part_of_lx_vector(char *out_e_or_o,char *in_lx,int bps,int par);
void take_e_part_of_lx_color(color *out_e,color *in_lx);
void take_e_part_of_lx_vector(char *out_e,char *in_lx,int bps);
void take_o_part_of_lx_color(color *out_o,color *in_lx);
void take_o_part_of_lx_vector(char *out_o,char *in_lx,int bps);
void unset_eo_geometry();
