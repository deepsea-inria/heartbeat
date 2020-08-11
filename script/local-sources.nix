let pkgs = import <nixpkgs> {}; in

let
  projectSrc = ../.;
  pbenchSrc = pkgs.fetchFromGitHub {
      owner  = "mikerainey";
      repo   = "pbench";
      rev    = "1c90259b594b6612bc6b9973564e89c297ad17b3";
      sha256 = "1440zavl3v74hcyg49h026vghhj1rv5lhfsb5rgfzmndfynzz7z0";
  };
  pbenchOcamlSrcs = import "${pbenchSrc}/nix/local-sources.nix";
in

let
  nixSrc = "${projectSrc}/script";
in

{

  nixSrc = nixSrc;
  
  pbenchSrc = pbenchSrc;

  pbenchOcamlSrcs = pbenchOcamlSrcs;
  
  benchOcamlSrc = "${projectSrc}/bench/";

  cmdlineSrc = ../../cmdline;

  cilkRtsSrc = ../../cilk-plus-rts-with-stats;

  chunkedseqSrc = ../../chunkedseq;

  sptlSrc = ../../sptl;

  pbbsIncludeSrc = ../../pbbs-include;

  cactusStackSrc = ../../cactus-stack;

  heartbeatSrc = ../.;
  
}
