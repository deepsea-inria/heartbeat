let pkgs = import <nixpkgs> {}; in

{

  cmdlineSrc = pkgs.fetchFromGitHub {
    owner  = "deepsea-inria";
    repo   = "cmdline";
    rev    = "c5f96b4aecb2019b5a690176195d37f7df3ed34b";
    sha256 = "1rz9bfdd5242gy3vq4n9vj2rcr5pwp0j4cjycpn3pm7rnwrrcjnh";
  };

  cilkRtsSrc = pkgs.fetchFromGitHub {
    owner  = "deepsea-inria";
    repo   = "cilk-plus-rts-with-stats";
    rev    = "ae723bd498d9dd09dcfadefceca092d7e28d4352";
    sha256 = "0bgsi9aqhjb80ds2c9mb1m3jzs6y3hwd2683f8gg301giifv17dn";
  };

  pbenchSrc = pkgs.fetchFromGitHub {
    owner  = "deepsea-inria";
    repo   = "pbench";
    rev    = "6d862442713cae93153c439a3e5e5867e95af958";
    sha256 = "1qdsl2amaqkckp2r2s4j5az8b0bpkbrgf8nc0zp575l1yjdzi2hh";
  };

  chunkedseqSrc = pkgs.fetchFromGitHub {
    owner  = "deepsea-inria";
    repo   = "chunkedseq";
    rev    = "d2925cf385fb43aff7eeb9c08cce88d321e5e02e";
    sha256 = "09qyv48vb2ispl3zrxmvbziwf6vzjh3la7vl115qgmkq67cxv78b";
  };

  sptlSrc = pkgs.fetchFromGitHub {
    owner  = "deepsea-inria";
    repo   = "sptl";
    rev    = "cb2bf943780efb63b4ab231f7bc5f3af631ed822";
    sha256 = "17rm2cz0qa2c3pi1ya3rv8143g3d7j0hgm9c9bd5q7chsf0abi3c";
  };

  pbbsIncludeSrc = pkgs.fetchFromGitHub {
    owner  = "deepsea-inria";
    repo   = "pbbs-include";
    rev    = "ef5ce72c4b4c26af78f7d91ea9b51336cd83a2e9";
    sha256 = "1g64j8gv6s9ggzhr2ky0y55s404cm0yrmdbhi5q4gfqaczbymyr4";
  };

  cactusStackSrc = pkgs.fetchFromGitHub {
    owner  = "deepsea-inria";
    repo   = "cactus-stack";
    rev    = "25d43a07a3adf7e0de84bc8f8551ba94d1d0a75b";
    sha256 = "129dxs4s922dj6hs3rqy9asdc91262z21mwkzsm4g4lcsxychffa";
  };
  
  heartbeatSrc = pkgs.fetchFromGitHub {
    owner  = "deepsea-inria";
    repo   = "heartbeat";
    rev    = "26de8102080f2c43914a50ff8ed930d3f4f33062";
    sha256 = "1ppfk1q71j6hq3mdc56519aa3vb0p0kfn10hblsrmaz09gqr488k";
  };

}
