# cmdline reverse complement function
function revcomp(){
  echo $1 | rev | tr ATGCatgcMRWSYKVHDBNmrwsykvhdbn TACGtacgKYWSRMBDHVNkywsrmbdhvn
}
