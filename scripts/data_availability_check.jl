# data availability check
varlist = (:t2m, :u10, :v10, :sp, :snr, :vpd_cf)

map(1950:2022) do yr
        map(varlist) do vn
                filelist = readdir("/Net/Groups/data_BGC/era5/e1/0d25_hourly/$vn/$yr/", join = true)
                filter!(!contains("back_extension"),filelist)
                filter!(endswith(".nc"),filelist)
                size(filelist)[1] < 12 ? println( "$vn $yr $(size(filelist))") : missing
        end
end

# v10 1951 (4,)
# v10 1952 (0,)