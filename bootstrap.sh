if [[ "$OSTYPE" == "linux-gnu"* ]]; then
        conda env create --prefix ./JUMPptm -f $PWD/jump_ptm_linux.yml
elif [[ "$OSTYPE" == "darwin"* ]]; then
        conda env create --prefix ./JUMPptm -f $PWD/jump_ptm_mac.yml
fi

$PWD/JUMPptm/bin/cpanm HTTP::Message~"<= 6.20" File::Copy File::Basename Scalar::Util LWP::UserAgent Set::Partition Sys::Hostname Spreadsheet::XLSX Statistics::Distributions
conda activate $PWD/JUMPptm
