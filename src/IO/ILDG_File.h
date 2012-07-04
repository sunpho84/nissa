#ifndef _ILDG_FILE_H
#define _ILDG_FILE_H
ILDG_File ILDG_File_open(char *path,int amode);
ILDG_File ILDG_File_open_for_read_only(char *path);
ILDG_Offset ILDG_File_get_position(ILDG_File &file);
ILDG_Offset ILDG_File_get_size(ILDG_File &file);
ILDG_header ILDG_File_get_next_record_header(ILDG_File &file);
bool ILDG_File_reached_EOF(ILDG_File &file);
bool get_MB_flag(ILDG_header &header);
bool get_ME_flag(ILDG_header &header);
int ILDG_File_search_record(ILDG_header &header,ILDG_File &file,const char *record_name);
void ILDG_File_close(ILDG_File &file);
void ILDG_File_read_all(void *data,ILDG_File &file,int nbytes_req);
void ILDG_File_read_ildg_data_all(void *data,ILDG_File &file,ILDG_header &header);
void ILDG_File_seek_to_next_eight_multiple(ILDG_File &file);
void ILDG_File_skip_nbytes(ILDG_File &file,ILDG_Offset nbytes);
void ILDG_File_skip_record(ILDG_File &file,ILDG_header header);
void ILDG_File_read_checksum(checksum check_read,ILDG_File &file);
#endif
