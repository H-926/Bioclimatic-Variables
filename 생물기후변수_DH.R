library(terra)
library(predicts)

# 경로
# 참고: sample_data/ 에는 참고용 샘플 1개씩만 포함되어 있으며, 실제 실행 시에는
# 전체 WorldClim 월별 자료(1990-2024) 및 dem.asc를 별도로 받아 아래 경로에 위치시켜야 함.
base   <- "WorldClim_10min_1990-2024"
dem    <- rast("dem.asc")
outdir <- "bcVariables9024"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# 월 평균 layer
mmean12 <- function(dirpath, prefix){
  months <- sprintf("%02d", 1:12)
  r <- lapply(months, function(m){
    f <- list.files(dirpath, pattern = paste0("^", prefix, "_\\d{4}-", m, "\\.tif$"),
                    full.names = TRUE)
    if(length(f)==0) stop("파일 없음: ", m)
    mean(rast(f), na.rm = TRUE)
  })
  x <- rast(r); names(x) <- months; x
}
tmin12 <- mmean12(file.path(base, "wc2.1_cruts4.09_10m_tmin_1990-2024"), "wc2.1_cruts4.09_10m_tmin")
tmax12 <- mmean12(file.path(base, "wc2.1_cruts4.09_10m_tmax_1990-2024"), "wc2.1_cruts4.09_10m_tmax")
prec12 <- mmean12(file.path(base, "wc2.1_cruts4.09_10m_prec_1990-2024"), "wc2.1_cruts4.09_10m_prec")

# WGS84 부여
if (is.na(crs(dem)))    crs(dem)    <- "+proj=longlat +datum=WGS84"
if (is.na(crs(tmin12))) crs(tmin12) <- "+proj=longlat +datum=WGS84"
if (is.na(crs(tmax12))) crs(tmax12) <- "+proj=longlat +datum=WGS84"
if (is.na(crs(prec12))) crs(prec12) <- "+proj=longlat +datum=WGS84"

# biovars 계산
bio <- predicts::bcvars(prec12, tmin12, tmax12)
names(bio) <- paste0("bio", 1:19)

# 투영 + 마스킹
bio_dem <- project(bio, dem, method = "bilinear")
bio_dem <- mask(bio_dem, dem)

# 저장
for(nm in names(bio_dem)){
  writeRaster(bio_dem[[nm]], file.path(outdir, paste0(nm, ".tif")),
              overwrite=TRUE, filetype="GTiff", gdal="COMPRESS=LZW", NAflag=-9999)

  asc_path <- file.path(outdir, paste0(nm, ".asc"))
  writeRaster(bio_dem[[nm]], asc_path,
              overwrite=TRUE, filetype="AAIGrid", NAflag=-9999,
              gdal=c("DECIMAL_PRECISION=6"))

  # 정수에 붙은 ".000..." 제거: -180.000000000000 → -180, -9999.000000 → -9999
  txt <- readLines(asc_path)
  txt <- gsub("(-?\\d+)\\.0+(?=\\s|$)", "\\1", txt, perl = TRUE)
  writeLines(txt, asc_path)

  cat("저장:", nm, "\n")
}
cat("완료 → ", normalizePath(outdir), "\n")
