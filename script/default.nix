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
  nbBuildCores ? 0,
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
    let docs = if buildDocs then [ pkgs.pandoc ] else []; in
    [ pbench sptl pbbs-include cmdline chunkedseq
      pkgs.makeWrapper pkgs.R pkgs.texlive.combined.scheme-small
      pkgs.ocaml gcc pkgs.wget
    ] ++ docs;
        
  configurePhase =
    let hwlocConfig =
      if useHwloc then ''
        USE_HWLOC=1
        USE_MANUAL_HWLOC_PATH=1
        MY_HWLOC_FLAGS=-I ${hwloc.dev}/include/
        MY_HWLOC_LIBS=-L ${hwloc.lib}/lib/ -lhwloc
      '' else "";
    in
    let settingsScript = pkgs.writeText "settings.sh" ''
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
      CILK_EXTRAS_PREFIX=-L ${cilk-plus-rts-with-stats}/lib -I ${cilk-plus-rts-with-stats}/include -ldl -DCILK_RUNTIME_WITH_STATS -DCUSTOM_CILK_PLUS_RUNTIME
      ${hwlocConfig}
    '';
    in
    let sptlConfigFile = pkgs.writeText "sptl_config.txt" "${sptl}/bin/"; in
    ''
    cp -r --no-preserve=mode ${pbench} pbench
    cp ${settingsScript} bench/settings.sh
    cp ${sptlConfigFile} bench/sptl_config.txt
    '';

  buildPhase =
    let getNbCoresScript = pkgs.writeScript "get-nb-cores.sh" ''
      #!/usr/bin/env bash
      ${sptl}/bin/get-nb-cores.sh
    '';
    in
    ''
    cp ${getNbCoresScript} bench/
    make -C bench bench.pbench
    '';  

  installPhase =
    let lu =
        if useLibunwind then
           ''--prefix LD_LIBRARY_PATH ":" ${libunwind}/lib''
        else "";
    in
    let hw =
        if useHwloc then
          ''--prefix LD_LIBRARY_PATH ":" ${hwloc.lib}/lib''
        else "";
    in
    let nmf = "-skip make";
    in
    let rf =
      if pathToResults != "" then
        "-path_to_results ${pathToResults}"
      else "";
    in
    let df =
      if pathToResults != "" then
        "-path_to_data ${pathToData}"
      else "";
    in
    let nbc =
      if nbBuildCores != 0 then
        "-nb_make_cores ${toString nbBuildCores}"
      else "";
    in
    let flags = "${nmf} ${rf} ${df} ${nbc}";
    in
    ''
    mkdir -p $out/bench/
    cp bench/bench.pbench bench/timeout.out $out/bench/
    wrapProgram $out/bench/bench.pbench --prefix PATH ":" ${pkgs.R}/bin \
       --prefix PATH ":" ${pkgs.texlive.combined.scheme-small}/bin \
       --prefix PATH ":" ${gcc}/bin \
       --prefix PATH ":" ${pkgs.wget}/bin \
       --prefix PATH ":" $out/bench \
       --prefix LD_LIBRARY_PATH ":" ${gcc}/lib \
       --prefix LD_LIBRARY_PATH ":" ${gcc}/lib64 \
       --prefix LD_LIBRARY_PATH ":" ${gperftools}/lib \
       --prefix LD_LIBRARY_PATH ":" ${cilk-plus-rts-with-stats}/lib \
       --set TCMALLOC_LARGE_ALLOC_REPORT_THRESHOLD 100000000000 \
       ${lu} \
       ${hw} \
       --add-flags "${flags}"
    pushd bench
    $out/bench/bench.pbench compare -only make
    popd
    cp bench/sptl_config.txt $out/bench/sptl_config.txt
    cp bench/*.heartbeat bench/*.cilk_elision bench/*.cilk $out/bench/
    '';

  meta = {
    description = "The C++ library for benchmarking the Heartbeat Scheduling algorithm for fine-grain, irregular parallelism.";
    license = "MIT";
    homepage = http://deepsea.inria.fr/heartbeat/;
  };
}
