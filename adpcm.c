// wav2adpcm transcoder
// - rlyeh, public domain

#include <stdio.h>
#include <stdlib.h>

#define JO_WAV_IMPLEMENTATION
#include "jo_wav.h"

#define DR_WAV_IMPLEMENTATION
#include "dr_wav.h"

int main(int argc, char **argv) {
    if( argc != 3 ) {
        return -fprintf(stderr, "%s infile.wav outfile.wav\n", argv[0]);
    }

    int rc = -1;
    char *infile = argv[1];
    char *outfile = argv[2];

    for( drwav w = {0}, *wav = &w; wav && drwav_init_file(wav, infile, NULL); wav = 0 ) {
        int channels = wav->channels;
        int frequency = wav->sampleRate;
        int numsamples = wav->totalPCMFrameCount;
        short *samples = realloc(0, sizeof(short) * numsamples * channels); // @leak
        if( samples ) drwav_read_pcm_frames_s16(wav, numsamples, samples);
        drwav_uninit(wav);

        // SLOW    Uses samples_per_block=2041, lookahead=3, greedy=0, quality=4
        // FAST    Uses samples_per_block=2041, lookahead=3, greedy=1, quality=4
        // FASTEST Uses samples_per_block=2041, lookahead=0, greedy=1, quality=0
        if( samples )
        rc = !jo_write_wav_adpcm_ex(outfile, channels, frequency, samples, numsamples,
                2041, // samples_per_block=2041: a good default. other common default is 505 (256 bytes per chunk). Must be multiple of 8 plus 1.
                3,    // lookahead=3: optimize this many steps in advance. 0 == faster, 1+ slower but better
                1,    // greedy=1: 0 .. 1 = faster, but less accurate
                4     // quality=4: 0 .. 4 (higher numbers are higher quality)
        );
    }

    return rc;
}
