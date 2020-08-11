let pkgs = import <nixpkgs> {}; in

let
  pbenchSrc = pkgs.fetchFromGitHub {
    owner  = "mikerainey";
    repo   = "pbench";
    rev    = "1c90259b594b6612bc6b9973564e89c297ad17b3";
    sha256 = "1440zavl3v74hcyg49h026vghhj1rv5lhfsb5rgfzmndfynzz7z0";
  };
  cmdlineSrc = pkgs.fetchFromGitHub {
    owner  = "deepsea-inria";
    repo   = "cmdline";
    rev    = "c5f96b4aecb2019b5a690176195d37f7df3ed34b";
    sha256 = "1rz9bfdd5242gy3vq4n9vj2rcr5pwp0j4cjycpn3pm7rnwrrcjnh";
  };
  chunkedseqSrc = pkgs.fetchFromGitHub {
    owner  = "deepsea-inria";
    repo   = "chunkedseq";
    rev    = "d2925cf385fb43aff7eeb9c08cce88d321e5e02e";
    sha256 = "09qyv48vb2ispl3zrxmvbziwf6vzjh3la7vl115qgmkq67cxv78b";
  };
  cilkRtsSrc = pkgs.fetchFromGitHub {
    owner  = "deepsea-inria";
    repo   = "cilk-plus-rts-with-stats";
    rev    = "d143c31554bc9c122d168ec22ed65e7941d4c91d";
    sha256 = "123bsrqcp6kq6xz2rn4bvj2nifflfci7rd9ij82fpi2x6xvvsmsb";
  };
  heartbeatSrc = ../.;
    # pkgs.fetchFromGitHub {
    #   owner  = "deepsea-inria";
    #   repo   = "heartbeat";
    #   rev    = "26de8102080f2c43914a50ff8ed930d3f4f33062";
    #   sha256 = "1ppfk1q71j6hq3mdc56519aa3vb0p0kfn10hblsrmaz09gqr488k";
  # };
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
  sptlSrc = pkgs.fetchFromGitHub {
    owner  = "deepsea-inria";
    repo   = "sptl";
    rev    = "911bc7af7c658020138a08d4923224332b08a27f";
    sha256 = "1w4nkcg1bcbkxr3i31ac9w29kly7ms5h8xkigss48n1hcizvwscz";
  };
in    

{

  cmdline = import "${cmdlineSrc}/script/default.nix" {};
  
  chunkedseq = import "${chunkedseqSrc}/script/default.nix" {};

  cilk-plus-rts-with-stats = import "${cilkRtsSrc}/default.nix" { };
  
  pbenchSrc = pbenchSrc;

  pbenchOcamlSrcs = import "${pbenchSrc}/nix/local-sources.nix";

  benchOcamlSrc = "${heartbeatSrc}/bench/";

  sptl = import "${sptlSrc}/script/default.nix" { };

  pbbs-include = import "${pbbsIncludeSrc}/default.nix" { };

  cactus-stack = import "${cactusStackSrc}/script/default.nix" { };
    
  heartbeatSrc = heartbeatSrc;
  
}
