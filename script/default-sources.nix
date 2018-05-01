let pkgs = import <nixpkgs> {}; in

{

  cmdlineSrc = ../../cmdline;

  cilkRtsSrc = ../../cilk-plus-rts-with-stats;

  pbenchSrc = ../../pbench;

  chunkedseqSrc = ../../chunkedseq;

  sptlSrc = ../../sptl;

  pbbsIncludeSrc = ../../pbbs-include;

  cactusStackSrc = ../../cactus-stack;

  heartbeatSrc = ../.;
  
}
