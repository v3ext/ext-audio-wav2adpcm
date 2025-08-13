#if CODE
AUTORUN {
    recipe(
        // infiles, final ziplevel
        "*.wav",0,
        // conversion steps
        "ext/ext-audio-wav2adpcm/adpcm.EXE INPUT OUTPUT.adpcm\n"
    );
}
#endif
