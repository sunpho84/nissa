#pragma once

void crash(const char *mess,int err)
{
  fprintf(stderr,"%s",mess);
  fflush(stderr);
  exit(1);
}
