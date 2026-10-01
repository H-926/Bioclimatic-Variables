# WorldClim 10min monthly data(1990-2024) -> bio1~19 (.tif + .asc)

library(terra)
library(predicts)

# 경로
# 참고: sample_data/ 에는 참고용 샘플 1개씩만 포함되어 있으며, 실제 실행 시에는
# 전체 WorldClim 월별 자료(1990-2024) 및 dem.asc를 별도로 받아 아래 경로에 위치시켜야 함.
base     <- "WorldClim_10min_1990-2024"
dem_path <- "dem.asc"
outdir   <- "bcVariables9024"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

YEARS   <- 1990:2024
N_YEAR  <- length(YEARS)
DEC     <- 2
SIGDIG  <- 7

stopifnot(dir.exists(base), file.exists(dem_path))
dem <- rast(dem_path)

mmean12 <- function(dirpath, prefix){
  months <- sprintf("%02d", 1:12)
  r <- lapply(months, function(m){
    f <- list.files(dirpath, pattern = paste0("^", prefix, "_\\d{4}-", m, "\\.tif$"),
                    full.names = TRUE)
    if (length(f) != N_YEAR)
      stop(sprintf("%s %s월: 파일 %d개, %d개 기대", basename(dirpath), m, length(f), N_YEAR))
    mean(rast(f), na.rm = TRUE)
  })
  x <- rast(r); names(x) <- months; x
}
tmin12 <- mmean12(file.path(base, "wc2.1_cruts4.09_10m_tmin_1990-2024"), "wc2.1_cruts4.09_10m_tmin")
tmax12 <- mmean12(file.path(base, "wc2.1_cruts4.09_10m_tmax_1990-2024"), "wc2.1_cruts4.09_10m_tmax")
prec12 <- mmean12(file.path(base, "wc2.1_cruts4.09_10m_prec_1990-2024"), "wc2.1_cruts4.09_10m_prec")
message("월 평균 layer 생성 완료 (tmin/tmax/prec x 12)")

for (nm in c("dem", "tmin12", "tmax12", "prec12")) {
  x <- get(nm)
  if (identical(crs(x), "")) { crs(x) <- "EPSG:4326"; assign(nm, x); message("CRS 부여: ", nm) }
}

bio <- predicts::bcvars(prec12, tmin12, tmax12)
stopifnot(nlyr(bio) == 19)
names(bio) <- paste0("bio", 1:19)

if (!compareGeom(bio, dem, stopOnError = FALSE))
  stop("bio와 dem의 격자가 다름. 월별 자료를 dem 격자로 resample한 뒤 bcvars 계산.")
bio_dem <- mask(bio, dem)

mx <- max(abs(c(global(bio_dem, "min", na.rm = TRUE)[, 1],
                global(bio_dem, "max", na.rm = TRUE)[, 1])), na.rm = TRUE)
if (!is.finite(mx) || mx <= 0) stop("전 레이어 최대 절대값이 유효하지 않음: ", mx)
need <- floor(log10(mx)) + 1 + DEC
if (need > SIGDIG) stop(sprintf("SIGNIFICANT_DIGITS=%d 부족: 최대값 %.4g에 %d자리 필요", SIGDIG, mx, need))
message(sprintf("최대 절대값 %.4g -> 필요 유효자릿수 %d (SIGDIG=%d)", mx, need, SIGDIG))

report <- data.frame()
for (nm in names(bio_dem)) {
  lyr <- bio_dem[[nm]]

  tif_path <- file.path(outdir, paste0(nm, ".tif"))
  writeRaster(lyr, tif_path, overwrite = TRUE, filetype = "GTiff",
              gdal = "COMPRESS=LZW", NAflag = -9999)

  asc_path <- file.path(outdir, paste0(nm, ".asc"))
  writeRaster(round(lyr, DEC), asc_path, overwrite = TRUE, filetype = "AAIGrid",
              NAflag = -9999, gdal = sprintf("SIGNIFICANT_DIGITS=%d", SIGDIG))

  n_line <- length(readLines(asc_path, warn = FALSE))
  if (n_line != 6 + nrow(lyr))
    stop(sprintf("%s: .asc %d줄, %d줄 기대 — 잘림", nm, n_line, 6 + nrow(lyr)))

  n_tif <- global(!is.na(lyr), "sum")[1, 1]
  n_asc <- global(!is.na(rast(asc_path)), "sum")[1, 1]
  if (n_asc != n_tif)
    stop(sprintf("%s: 유효셀 .asc %d vs .tif %d", nm, n_asc, n_tif))

  report <- rbind(report, data.frame(
    layer = nm, valid = n_tif, asc_MB = round(file.size(asc_path) / 1e6, 2),
    tif_MB = round(file.size(tif_path) / 1e6, 2)))
  cat("저장:", nm, "\n")
}

print(report, row.names = FALSE)
message(sprintf("합계 .asc %.0f MB / .tif %.0f MB", sum(report$asc_MB), sum(report$tif_MB)))
message("완료 -> ", normalizePath(outdir))
