#!/bin/sh

gcc calc_ground_state.c -lm -o calc_ground_state.out -Wall

./calc_ground_state.out > dat

rm calc_ground_state.out
