{ pkgs   ? import <nixpkgs> {},
  stdenv ? pkgs.stdenv,
  sources ? import ./default-sources.nix,
  gperftools ? pkgs.gperftools,
  useHwloc ? false,
  hwloc ? pkgs.hwloc,
  libunwind ? pkgs.libunwind,
  useLibunwind ? false,
  gcc ? pkgs.gcc,
  pathToResults ? "",
  pathToData ? "",
  buildDocs ? false
}:

let

  callPackage = pkgs.lib.callPackageWith (pkgs // sources // self);

  self = {

    hwloc = hwloc;
    useHwloc = useHwloc;

    libunwind = libunwind;
    useLibunwind = useLibunwind;

    gperftools = gperftools;

    gcc = gcc;

    buildDocs = buildDocs;

    pbench = callPackage "${sources.pbenchSrc}/script/default.nix" { };
    cmdline = callPackage "${sources.cmdlineSrc}/script/default.nix" { };
    cilk-plus-rts-with-stats = callPackage "${sources.cilkRtsSrc}/script/default.nix" { };
    chunkedseq = callPackage "${sources.chunkedseqSrc}/script/default.nix" { };
    pbbs-include = callPackage "${sources.pbbsIncludeSrc}/default.nix" { };
    sptl = callPackage "${sources.sptlSrc}/script/default.nix" { };
    cactus-stack = callPackage "${sources.cactusStackSrc}/script/default.nix" { };
    heartbeatSrc = sources.heartbeatSrc;
    
  };

in

with self;

stdenv.mkDerivation rec {
  name = "heartbeat-${version}";
  version = "v1.0-pldi";

  src = heartbeatSrc;

  buildInputs =
    let docs =
      if buildDocs then [
        pkgs.pandoc
        pkgs.texlive.combined.scheme-full
      ] else
        [];
    in
    [ pbench sptl pbbs-include cmdline chunkedseq ] ++ docs;
        
  buildPhase =
    let hwlocConfig =
      if useHwloc then ''
        USE_HWLOC=1
        USE_MANUAL_HWLOC_PATH=1
        MY_HWLOC_FLAGS=-I ${hwloc.dev}/include/
        MY_HWLOC_LIBS=-L ${hwloc.lib}/lib/ -lhwloc
      '' else "";
    in
    ''
    cat >> settings.sh <<__EOT__
    CACTUS_PATH=${cactus-stack}/include/
    PBENCH_PATH=../pbench/
    CMDLINE_PATH=${cmdline}/include/
    CHUNKEDSEQ_PATH=${chunkedseq}/include/
    SPTL_PATH=${sptl}/include/
    PBBS_INCLUDE_PATH=${pbbs-include}/include/
    ENCORE_INCLUDE_PATH=$out/include/
    USE_32_BIT_WORD_SIZE=1
    USE_CILK=1
    CUSTOM_MALLOC_PREFIX=-ltcmalloc -L${gperftools}/lib
    CILK_EXTRAS_PREFIX=-L ${cilk-plus-rts-with-stats}/lib -I  ${cilk-plus-rts-with-stats}/include -ldl -DCILK_RUNTIME_WITH_STATS
    __EOT__
    cat >> settings.sh <<__EOT__
    ${hwlocConfig}
    __EOT__
    '';

  installPhase = ''
    mkdir -p $out/bench/
    cp settings.sh bench/Makefile bench/bench.ml bench/*.cpp bench/*.hpp $out/bench/
    mkdir -p $out/include/
    cp include/*.hpp $out/include/
    '';

  meta = {
    description = "The C++ library for benchmarking the Heartbeat Scheduling algorithm for fine-grain, irregular parallelism.";
    license = "MIT";
    homepage = http://deepsea.inria.fr/heartbeat/;
  };
}
