#!/bin/bash
#SBATCH --time=0-6:00:00
#SBATCH --account=def-jbl
#SBATCH --mem=64G
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --gres=gpu:1
#SBATCH --job-name=loocv_subset_bam
#SBATCH --output=logs/%x-%A_fold%a.out
#SBATCH --error=logs/%x-%A_fold%a.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=sophie.bernard-doucet@umontreal.ca
#SBATCH --array=0-3

#""" 10x genomics, command from email """