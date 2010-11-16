# @ shell = /bin/bash
# @ job_name = job_pirl
# @ error  = $(job_name).$(jobid).err
# @ output = $(job_name).$(jobid).out
# @ environment = COPY_ALL
# @ notification = always
# @ notify_user = Francesco.Sanfilippo@roma1.infn.it
# @ wall_clock_limit = 03:00:00
# @ job_type         = BLUEGENE
# @ bg_size          = 512
# @ bg_connection    = TORUS
# @ queue
#
