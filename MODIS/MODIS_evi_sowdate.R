library(terra)
library (zoo)
library (signal)


# allrasters1<-stack(allrasters)
get_sgolayfilt<-function(x){
  if (!all(is.na(x))) {
    v <- as.vector(x)
    z <- zoo::na.spline(v)
    s1.ts2 <- ts(z,start=1, end=length(x))
    # as.numeric(s1.ts2[1])
    x <- sgolayfilt(s1.ts2)
  }
  return(x)
}
# fun <- function(x) {
#   mean(x)
# }

# extract(firstsgol, 1331)
# sgol_df <- as.data.frame(firstsgol)
# x <- as.numeric(as.vector(sgol_df[1331, ]))

get_sowday_from_sgolay <- function(x, first_image_doy, last_image_doy) {
  evi <- x * 0.0001
  max_evi <- max(evi)
  sowday <- -9999
  # if (max_evi < 0.5) {
  #   sowday <- 0
  # } else {
  # start_doy <- as.numeric(strftime(as.Date("2019-04-01"),format = "%j"))
  # end_doy <- as.numeric(strftime(as.Date("2019-12-31"),format = "%j"))
  tryCatch(
    {
      days <- seq(first_image_doy, last_image_doy, length.out = length(evi))

      # determine index of max value of evi
      max_idx <- which.max(evi)[1]
      day_max <- days[max_idx]

      # shift EVI so that it contains a minimum value of 0 (prior to the)
      evi_0shift <- evi - min(evi[1:max_idx])

      # get sow date: day where evi is equal to 10% of max evi
      evi_10pct_max_shift <- 0.1 * max(evi_0shift)
      evi_low_idx <- which(evi_0shift < evi_10pct_max_shift)
      min_idx <- max(evi_low_idx[evi_low_idx < max_idx])

      x1 <- days[min_idx]
      x2 <- days[min_idx + 1]
      y1 <- evi_0shift[min_idx]
      y2 <- evi_0shift[min_idx + 1]
      day_10pct <- (evi_10pct_max_shift - y1) * (x2 - x1) / (y2 - y1) + x1


      # get crop maturity date: day where evi is equal to 10% less than max evi
      evi_90pct_max_shift <- 0.9 * max(evi_0shift)
      evi_high_idx <- which(evi_0shift < evi_90pct_max_shift)
      post_maturity_idx <- min(evi_high_idx[evi_high_idx > max_idx])

      x1 <- days[post_maturity_idx - 1]
      x2 <- days[post_maturity_idx]
      y1 <- evi_0shift[post_maturity_idx - 1]
      y2 <- evi_0shift[post_maturity_idx]
      day_90pct <- (evi_90pct_max_shift - y1) * (x2 - x1) / (y2 - y1) + x1

      # # plot key dates and evi_shift
      # plot(days, evi_0shift)
      # lines(days, rep(evi_10pct_max_shift, length(evi_0shift)), col = 'red')
      # lines(days, rep(evi_90pct_max_shift, length(evi_0shift)), col = 'red')
      # points(c(day_10pct, day_max, day_90pct),
      #        c(evi_10pct_max_shift, max(evi_0shift), evi_90pct_max_shift),
      #        col = 'red')

      # get growing season length
      day_sow <- day_10pct - 15
      day_maturity <- day_90pct
      growing_season_days <- day_90pct - day_sow

      # only count as rice if the growing season is between 112 and 152 days
      # if (growing_season_days < 112 | growing_season_days > 152) {
      #   sowday <- -growing_season_days
      # } else {
      #   sowday <- day_sow
      # }

      sowday <- day_sow
    },
    error = function(cond) {
      sowday <- -2
      growing_season_days <- -2
    },
    warning = function(cond) {
      sowday <- -3
      growing_season_days <- -3
    },
    finally = {
      # do nothing
      if (!exists('growing_season_days')) {
        growing_season_days <- as.numeric(NA)
      }
    }
  )
  # }
  if (sowday == -9999) {
    sowday <- as.numeric(NA)
  }
  sow_calcs <- c(sowday, max_evi, growing_season_days)
  names(sow_calcs) <- c("sowday", "max_evi", "growing_season_days")
  return(sow_calcs)
}


get_year_sowday <- function(rasters_year, rastdates, out_path) {

  if (file.exists(out_path)) {
    cat('The file already exists. Reading from ', out_path, '...')
    sowdays <- rast(out_path)
    # sowdays$sowday[sowdays$sowday==-9999] <- NA
    # writeRaster(sowdays, gsub('\\.tif','_2.tif',out_path), overwrite = TRUE)
    cat('done.\n')
  } else {
    first_image_doy <- min(as.numeric(strftime(rastdates,"%j")))
    last_image_doy <- max(as.numeric(strftime(rastdates,"%j")))

    cat('Running Savitsky-golay filter for all pixels...')
    firstsgol<-terra::app(rasters_year, get_sgolayfilt)
    cat('done.\n')

    cat('Getting sow date and growing season for all pixels...')
    sowdays<-terra::app(firstsgol, get_sowday_from_sgolay,
                        first_image_doy = first_image_doy, last_image_doy = last_image_doy)
    cat('done.\n')


    # # filter for max EVI > 0.5
    # sowdays[sowdays$max_evi < 0.5] <- NA
    # # remove NA values
    # sowdays$sowday[sowdays$sowday==-9999] <- NA

    cat("Writing to file ", out_path)
    writeRaster(sowdays, out_path)
    cat('done.\n')
  }

  return(sowdays)
}



library(sf)
library(tidyverse)
subdists_india <- read_sf('/Users/gopal/Projects/Data/India/admin_boundaries/SUBDISTRICT_11.shp')
subdists_punjab <- subdists_india %>% dplyr::filter(STATE_UT == 'Punjab')

subdists_punjab_vect <- vect(subdists_punjab) #, crs = crs(year_sowday, proj = TRUE)
subdists_punjab_vect <- project(subdists_punjab_vect, year_sowday)

punjab_evi_path <- '/Users/gopal/Google Drive/_Research/Research projects/G_D/GEE_Div/GEE_DIV_MODIS_EVI_punjab/GEE_DIV_MODIS_EVI_punjab_evi'

rastlist <- list.files(path = punjab_evi_path, pattern='.tif$',
                       all.files=TRUE, full.names=T)
rastdates <- as.Date(gsub(".*?([0-9][0-9][0-9][0-9]_[0-9][0-9]_[0-9][0-9])\\.tif$","\\1",rastlist),"%Y_%m_%d")
rast_years <-as.numeric(strftime(rastdates, "%Y"))


out_dir <- file.path(punjab_evi_path, "sowdates")
if (!file.exists(out_dir)) dir.create(out_dir)


median_na_rm <- function(x) median(x, na.rm = TRUE)


set.seed(101)
px_samples <- sample(1:(nrow(year_sowday)*ncol(year_sowday)), size = 40, replace = FALSE)

years <- 2019:2022
# year_i <- 2019
for (year_i in years) {

  cat('getting year',year_i)
  rastdates_in_year <- rastdates[rast_years == year_i]
  rastlist_in_year <- rastlist[rast_years == year_i]
  rasters_year<-rast(rastlist_in_year)

  cat('removing zeros from EVI timeseries...')
  for (i in 1:(dim(rasters_year)[3])) {
    rasters_year[[i]][rasters_year[[i]] == 0] <- NA
  }
  cat('done.\n')

  out_path <- file.path(out_dir, paste0("sowdate_",year_i,".tif"))

  year_sowday <- get_year_sowday(rasters_year, rastdates_in_year, out_path)

  year_sowday[year_sowday$max_evi < 0.5] <- NA

  # png(file.path(out_dir, paste0('plot_sowdate_',year_i,".png")), width = 1200, height = 1000)
  # plot(year_sowday)
  # dev.off()

  subdist_sowday <- terra::extract(year_sowday[['sowday']], subdists_punjab_vect, fun = median_na_rm)$sowday
  if (year_i == years[1]) {
    subdist_sowday_df <- data.frame(SUB_ID = subdists_punjab$SUB_ID, sowday = subdist_sowday, year = year_i)
  } else {
    subdist_sowday_df <- subdist_sowday_df %>%
      dplyr::bind_rows(data.frame(SUB_ID = subdists_punjab$SUB_ID, sowday = subdist_sowday, year = year_i))
  }

  # plot intraannual timeseries
  ts_sowdate <- tibble(id = px_samples) %>%
    bind_cols(year_sowday[px_samples], year = year_i) %>%
    mutate(date = as.Date(paste0(year, "-", sowday), "%Y-%j"))
  ts_evi <- tibble(id = px_samples) %>% bind_cols(rasters_year[px_samples]) %>%
    setNames(c("id", as.character(rastdates_in_year))) %>%
    pivot_longer(starts_with("2"), names_to = "date", values_to = "evi") %>%
    mutate(date = as.Date(date))

  p_ts_growing_season <- ggplot(ts_evi) +
    geom_line(aes(date, evi)) +
    geom_segment(data = ts_sowdate, aes(date, 5000, xend = date + growing_season_days, yend = 5000), color = 'red') +
    facet_wrap(~id)
  ggsave(paste0('growing_season_ts_',year_i,'.png'), p_ts_growing_season, path = out_dir, width = 12, height = 12)
}

write_csv(subdist_sowday_df, file.path(out_dir,"subsdist_median_sowday_2019_2022.csv"))

subdists_sowday_sf <- subdists_punjab %>%
  full_join(subdist_sowday_df, by= "SUB_ID")
subdists_sowday_sf$sowdate <- as.Date(paste0("0000-", subdists_sowday_sf$sowday), "%Y-%j")

p_subdist_sowday_maps <- ggplot() + geom_sf(data = subdists_sowday_sf, aes(fill = sowdate)) +
  facet_wrap(~year)
ggsave("subdist_sowday_maps.png", p_subdist_sowday_maps, path = out_dir, width = 6, height = 6)



# orig_df <- as.data.frame(allrasters)
# sgol_df <- as.data.frame(firstsgol)
#
# set.seed(100)
# px_idx <- sample(1:nrow(orig_df), 20, replace = FALSE)
# ts_df <- tibble(px_idx = px_idx) %>%
#   bind_cols(orig_df[px_idx, ] %>% setNames(paste0("orig_", 1:ncol(orig_df)))) %>%
#   bind_cols(sgol_df[px_idx, ] %>% setNames(paste0("sgolay_", 1:ncol(sgol_df)))) %>%
#   pivot_longer(matches("(orig.*)|(sgolay.*)"), names_to = "obs_id", values_to="EVI") %>%
#   separate_wider_delim(cols = "obs_id", names = c("type","obs_num"), delim = "_") %>%
#   mutate(obs_num = as.numeric(obs_num))
#
#
# ggplot(ts_df, aes(obs_num, EVI, color = type)) +
#   geom_point() + geom_line() +
#   facet_wrap(~px_idx)
