#!/bin/bash

for i in $(cat config/timestamp_ordered_list)
do
    touch $i
done