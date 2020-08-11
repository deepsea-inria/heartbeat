# to build w/ clang:
# $ nix-shell --arg clang '(import <nixpkgs> {}).clang'

{ pkgs   ? import <nixpkgs> {},
  stdenv ? pkgs.stdenv,
  sources ? import ../script/default-sources.nix,
  cmdline ? sources.cmdline,
  cilk-plus-rts-with-stats ? sources.cilk-plus-rts-with-stats,
  chunkedseq ? sources.chunkedseq,
  cactus-stack ? sources.cactus-stack,
  sptl ? sources.sptl,
  pbbs-include ? sources.pbbs-include,
  gperftools ? pkgs.gperftools,
  gcc ? pkgs.gcc7,
  clang ? null,
  hwloc ? pkgs.hwloc, # use hwloc, unless this parameter equals null
  pviewSrc ? pkgs.fetchFromGitHub {
    owner  = "deepsea-inria";
    repo   = "pview";
    rev    = "78d432b80cc1ea2767e1172d56e473f484db7f51";
    sha256 = "1hd57237xrdczc6x2gxpf304iv7xmd5dlsvqdlsi2kzvkzysjaqn";
  }
}:

with pkgs; {
  qpidEnv = stdenvNoCC.mkDerivation {
    name = "cpp";
    
    buildInputs =
      [ gcc gperftools which cilk-plus-rts-with-stats]
      ++ (if hwloc == null then [] else [ hwloc ])
      ++ (if clang != null then [ clang ] else []);
    
    shellHook =
      let clangConfig =
          if clang != null then ''
          export CXX="${clang}/bin/clang++ -frtti -fno-exceptions"
          export LD_FLAGS=-latomic $LD_FLAGS
          '' else ''
            export CXX="${gcc}/bin/c++"
          '';
      in
      let hwlocFlgs =
            if hwloc == null then
              ""
            else
              ''export USE_HWLOC=1
                export HWLOC_FLAGS="-I ${hwloc.dev}/include/"
                export HWLOC_LIBS="-L ${hwloc.lib}/lib/ -lhwloc"
              '';
      in
      let pview = import "${pviewSrc}/default.nix" {}; in
      ''
      ${clangConfig}
      export CMDLINE_PATH="${cmdline}/include"
      export CILK_EXTRAS_PREFIX="-L ${cilk-plus-rts-with-stats}/lib -I ${cilk-plus-rts-with-stats}/include -ldl -DCILK_RUNTIME_WITH_STATS -DCUSTOM_CILK_PLUS_RUNTIME"
      export CHUNKEDSEQ_PATH="${chunkedseq}/include/"
      export SPTL_PATH="${sptl}/include/"
      export PBBS_INCLUDE_PATH="${pbbs-include}/include/"
      export HEARTBEAT_INCLUDE_PATH="../include/"
      export CACTUS_PATH="${cactus-stack}/include/"
      export USE_CILK=1
      export CUSTOM_MALLOC_PREFIX="-ltcmalloc -L ${gperftools}/lib"
      ${hwlocFlgs}
      export PATH=${gcc}/bin/:${pview}/bin:$PATH
      '';
  };
}
