#!/bin/bash

ls -rt $(cat config/keep_orderd_timestamp_list) > config/timestamp_ordered_list