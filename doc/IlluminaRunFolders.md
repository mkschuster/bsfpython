# Illumina Run Folders

## Data/Intensities/BaseCalls

 - One *.bcl file per lane (L001 - L008), cycle (C1.1) and tile
 - One *.filter file per lane and tile

### HiSeq 3000/4000

On a HiSeq 3000 with 310 cycles this amounts to 310 * 112 + 122 = 34,832 files per lane or
8 * 310 * 112 + 8 * 112 = 278,656 files in total in 8 * 310 + 8 = 2,488 directories.

## Data/Intensities/BaseCalls/Matrix

This directory exists only on the HiSeq 3000 and 4000 platforms.

 - One Matrix file per lane, cycle and tile

### HiSeq 3000/4000

On a HiSeq 3000 with 310 cycles this amounts to 310 * 112 = 34,720 files per lane or
8 * 310 * 112 = 277,760 files in total in 8 * 310 + 8 = 2,488 directories.

## Thumbnail_Images

 - One image file per lane, cycle, tile and channel (ACGT)
 - Two files (*.jpg and *.jpg.zprof) per surface, swath and channel

### HiSeq 2000

On a HiSeq 2000 with 209 cycles this amounts to 8 * 209 * 96 * 4 + 8 * 209 * 2 * 3 * 4 * 2 =
642,048 + 80,256 = 722,304 files in 8 * 209 + 8 = 1,680 directories.

## Notes

The whole HiSeq 3000 folder is 521 GiB and could be archived as a single file.
