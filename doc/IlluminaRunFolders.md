# Illumina Run Folders

## Software Versions

| RTA       | Application                     | Instrument      |
|:----------|:--------------------------------|:----------------|
| 1.12.4    | HiSeq Control Software 1.4.5    |                 |
| 1.12.4.2  | HiSeq Control Software 1.4.8    |                 |
| 1.13.48   | HiSeq Control Software 1.5.15.1 |                 |
| 1.17.21.3 | HiSeq Control Software 2.0.12.0 |                 |
| 1.18.54   | MiSeq Control Software 2.5.0.5  | MiSeq           |
| 1.18.61   | HiSeq Control Software 2.2.38   |                 |
| 1.18.64   | HiSeq Control Software 2.2.58   |                 |
| 2.5.2     | HiSeq Control Software 3.3.20   | HiSeq 3000/4000 |
| 2.7.3     | HiSeq Control Software 3.3.52   | HiSeq 3000/4000 |

- `IRF/Data/ImageSize.dat` not in RTA 2.5.2 and 2.7.3
- `IRF/Data/Intensities/config.xml` not in RTA 2.5.2 and 2.7.3
- `IRF/Data/Intensities/BaseCalls/config.xml` not in RTA 2.5.2 and 2.7.3
- `IRF/Data/Intensities/BaseCalls/L00{lane}/`
- `IRF/Data/Intensities/BaseCalls/L00{lane}/C{cycle}.1/`
- `IRF/Data/Intensities/BaseCalls/L00{lane}/C{cycle}.1/s_{lane}_{tile}.bcl(.gz)`
- `IRF/Data/Intensities/BaseCalls/L00{lane}/C{cycle}.1/s_{lane}_{tile}.stats` not in RTA 2.5.2 and 2.7.3
- `IRF/Data/Intensities/BaseCalls/L00{lane}/s_{lane}_{tile}.control` not in RTA 2.5.2 and 2.7.3
- `IRF/Data/Intensities/BaseCalls/L00{lane}/s_{lane}_{tile}.filter`
- `IRF/Data/Intensities/BaseCalls/Matrix/` not in RTA 2.7.3
- `IRF/Data/Intensities/BaseCalls/Matrix/L00{lane}` only in RTA 2.5.2
- `IRF/Data/Intensities/BaseCalls/Matrix/L00{lane}/C{cycle}.1`
- `IRF/Data/Intensities/BaseCalls/Matrix/L00{lane}/C{cycle}.1/s_{lane}_{tile}_matrix.txt`
- `IRF/Data/Intensities/BaseCalls/Matrix/s_{read}_matrix.txt` not in RTA 2.5.2
- `IRF/Data/Intensities/BaseCalls/Matrix/s_{lane}_{read}_matrix.txt` not in RTA 2.5.2
- `IRF/Data/Intensities/BaseCalls/Matrix/s_{lane}_{read}_{tile}_matrix.txt` not in RTA 2.5.2
- `IRF/Data/Intensities/BaseCalls/Phasing/` not in RTA 2.7.3
- `IRF/Data/Intensities/BaseCalls/Phasing/s_{read}_phasing.txt`
- `IRF/Data/Intensities/BaseCalls/Phasing/s_{lane}_{read}_phasing.txt`
- `IRF/Data/Intensities/BaseCalls/Phasing/s_{lane}_{read}_{tile}_cycle.txt` only in RTA 1.12.4, 1.12.4.2, 1.13.48 and
  1.17.21.3
- `IRF/Data/Intensities/BaseCalls/Phasing/s_{lane}_{read}_{tile}_phasing.txt`
- `IRF/Data/Intensities/BaseCalls/Phasing/s_{lane}_{cycle}_phasing.xml` not in RTA 2.5.2
- `IRF/Data/Intensities/BaseCalls/Phasing/EmpiricalPhasingCorrection_{lane}_{read}_{tile}.txt` not in 1.12.4, 1.12.4.2,
  1.13.48 and 1.17.21.3
- `IRF/Data/Intensities/BaseCalls/SampleSheet.csv` only in RTA 1.18.54
- `IRF/Data/Intensities/L00{lane}/` not in RTA 2.5.2 and 2.7.3
- `IRF/Data/Intensities/L00{lane}/s_{lane}_{tile}.locs` only in RTA 1.18.54
- `IRF/Data/Intensities/L00{lane}/s_{lane}_{tile}.clocs`
- `IRF/Data/Intensities/L00{lane}/C{cycle}.1/` not in RTA 1.18.54 and 1.18.64
- `IRF/Data/Intensities/L00{lane}/C{cycle}.1/s_{lane}_{tile}.cif`
- `IRF/Data/Intensities/L00{lane}/C{cycle}.1/s_{lane}_{tile}.errorMap` only in RTA 1.12.4, 1.12.4.2 and 1.13.48
- `IRF/Data/Intensities/L00{lane}/C{cycle}.1/s_{lane}_{tile}_T.FWHMMap` only in RTA 1.12.4, 1.12.4.2 and 1.13.48
- `IRF/Data/Intensities/Offsets/` not in RTA 2.5.2 and 2.7.3
- `IRF/Data/Intensities/Offsets/offsets.txt`
- `IRF/Data/Intensities/Offsets/SubTileOffsets.txt`
- `IRF/Data/Intensities/RTAConfiguration.xml` not in RTA 2.5.2 and 2.7.3
- `IRF/Data/Intensities/s_{lane}_{tile}_pos.txt` only in RTA 1.12.4, 1.12.4.2 and 1.13.48
- `IRF/Data/Intensities/s.locs` only in RTA 2.5.2 and 2.7.3
- `IRF/Data/reports/` only in RTA 1.12.4, 1.12.4.2 and 1.13.48
- `IRF/Data/RTALogs/` not in RTA 2.5.2 and 2.7.3
- `IRF/Data/Status_Files/` only in RTA 1.12.4, 1.12.4.2 and 1.13.48
- `IRF/Data/Status.htm` only in RTA 1.12.4, 1.12.4.2 and 1.13.48
- `IRF/Data/TileStatus/` only in RTA 1.18.54
- `IRF/Data/TileStatus/TileStatusL{lane}T{tile}.bin`
- `IRF/Data/TileStatus/TileStatusL{lane}T{tile}.tpl`
- `IRF/InterOp/`
- `IRF/InterOp/ColorMatrixMetricsOut.bin` only in RTA 2.7.3
- `IRF/InterOp/ControlMetricsOut.bin` not in RTA 2.5.2 and 2.7.3
- `IRF/InterOp/CorrectedIntMetricsOut.bin`
- `IRF/InterOp/EmpiricalPhasingMetricsOut.bin` only in RTA 2.5.2 and 2.7.3
- `IRF/InterOp/ErrorMetricsOut.bin`
- `IRF/InterOp/ExtractionMetricsOut.bin`
- `IRF/InterOp/EventMetricsOut.bin` only in RTA 2.5.2 and 2.7.3
- `IRF/InterOp/FWHMGridMetricsOut.bin` only in RTA 2.7.3
- `IRF/InterOp/ImageMetricsOut.bin` not in RTA 1.18.54, 2.5.2 and 2.7.3
- `IRF/InterOp/ImageMetricsOut.bin` only in RTA 2.7.3
- `IRF/InterOp/IndexMetricsOut.bin` only in RTA 1.18.54
- `IRF/InterOp/PFGridMetricsOut.bin` only in RTA 2.5.2 and 2.7.3
- `IRF/InterOp/QMetricsOut.bin`
- `IRF/InterOp/RegistrationMetricsOut.bin` only in RTA 2.5.2 and 2.7.3
- `IRF/InterOp/StaticRunMetricsOut.bin` only in RTA 2.7.3
- `IRF/InterOp/TileMetricsOut.bin`
- `IRF/PeriodicSaveRates/`
- `IRF/PeriodicSaveRates/Save All Thumbnails.xml`
- `IRF/Recipe/`
- `IRF/Recipe/{FlowCellIdentifier}.xml` not in RTA 1.18.54
- `IRF/Recipe/{ReagentKitRFIDTag/SerialNumber}.xml` only in RTA 1.18.54
- `IRF/Recipe/RunState.xml`
- `IRF/Recipe/{FlowCellIdentifier}_RunState.xml` not in RTA 1.18.54, 2.5.2 and 2.7.3
- `IRF/Thumbnail_Images/`
- `IRF/Thumbnail_Images/L00{lane}/`
- `IRF/Thumbnail_Images/L00{lane}/C{cycle}.1/`

## Data/Intensities/BaseCalls

- One *.bcl file per lane (L001 - L008), cycle (C1.1) and tile
- One *.filter file per lane and tile

### HiSeq 3000/4000

On a HiSeq 3000 with 310 cycles this amounts to 310 * 112 + 112 = 34,832 files per lane or
8 * 310 * 112 + 8 * 112 = 278,656 files in total in 8 * 310 + 8 = 2,488 directories.

## Data/Intensities/BaseCalls/Matrix

- One Matrix file per lane, cycle and tile

### HiSeq 3000/4000

On a HiSeq 3000 with 310 cycles this amounts to 310 * 112 = 34,720 files per lane or
8 * 310 * 112 = 277,760 files in total in 8 * 310 + 8 = 2,488 directories.

### HiSeq 2000

This directory exists only on the HiSeq 3000 and 4000 platforms.

## Data/Intensities/BaseCalls/Phasing

- One file per read (s_1_phasing.txt, s_2_phasing.txt, ...)
- One file per lane and read (s_1_1_phasing.txt, ..., s_8_3_phasing.txt)
- One file per surface, swath and tile (s_1_1_1101_phasing.txt, ..., s_8_3_2316_phasing.txt)

### HiSeq 3000/4000

This directory does not exist on this platform.

## Thumbnail_Images

- One image file per lane, cycle, tile and channel (ACGT)
- Two files (*.jpg and *.jpg.zprof) per surface, swath and channel

### HiSeq 2000

On a HiSeq 2000 with 209 cycles this amounts to 8 * 209 * 96 * 4 + 8 * 209 * 2 * 3 * 4 * 2 =
642,048 + 80,256 = 722,304 files in 8 * 209 + 8 = 1,680 directories.

## Notes

The whole HiSeq 3000 folder is 521 GiB and could be archived as a single file.
