//< @r-lyeh: modified jo_write_wav() to have float support (use `-32` bits)

// Written by Jon Olick 
//
// This is free and unencumbered software released into the public domain.
//
// Anyone is free to copy, modify, publish, use, compile, sell, or
// distribute this software, either in source code form or as a compiled
// binary, for any purpose, commercial or non-commercial, and by any
// means.
//
// In jurisdictions that recognize copyright laws, the author or authors
// of this software dedicate any and all copyright interest in the
// software to the public domain. We make this dedication for the benefit
// of the public at large and to the detriment of our heirs and
// successors. We intend this dedication to be an overt act of
// relinquishment in perpetuity of all present and future rights to this
// software under copyright law.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
// EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
// MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
// IN NO EVENT SHALL THE AUTHORS BE LIABLE FOR ANY CLAIM, DAMAGES OR
// OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE,
// ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
// OTHER DEALINGS IN THE SOFTWARE.
//
// For more information, please refer to <http://unlicense.org/>
//
// Inspiration from ADPCM-XQ 
//
// For implementation: 
//  #define JO_WAV_IMPLEMENTATION
//

#ifdef __cplusplus
extern "C" {
#endif

#ifndef JO_INCLUDE_WAV_H
#define JO_INCLUDE_WAV_H

int jo_write_wav(const char *filename, short num_channels, short bits_per_sample, int sample_rate_hz, const void *samples, int size_bytes);

int jo_write_wav_adpcm_ex(const char *filename, short num_channels, int sample_rate_hz, const short *samples, int num_samples, 
        int samples_per_block/*=2041*/,                     // a good default. other common default is 505 (256 bytes per chunk). Must be multiple of 8 plus 1.
        int lookahead/*=3*/,                                // optimize this many steps in advance. 0 == faster, 1+ slower but better
        int greedy/*=1*/,                                   // 0 .. 1 = faster, but less accurate
        int quality/*=4*/                                   // 0 .. 4 (higher numbers are higher quality)
        );

// Uses samples_per_block=2041, lookahead=0, greedy=1, quality=0
int jo_write_wav_adpcm_fastest(const char *filename, short num_channels, int sample_rate_hz, const short *samples, int num_samples);

// Uses samples_per_block=2041, lookahead=3, greedy=1, quality=4
int jo_write_wav_adpcm_fast(const char *filename, short num_channels, int sample_rate_hz, const short *samples, int num_samples);

// Uses samples_per_block=2041, lookahead=3, greedy=0, quality=4
int jo_write_wav_adpcm_slow(const char *filename, short num_channels, int sample_rate_hz, const short *samples, int num_samples);

#endif // JO_INCLUDE_WAV_H

//#define JO_WAV_IMPLEMENTATION
#ifdef JO_WAV_IMPLEMENTATION

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int jo_write_wav(const char *filename, short num_channels, short bits_per_sample, int sample_rate_hz, const void *samples, int size_bytes) {
    FILE *fp = fopen(filename, "wb");
    if(!fp) {
        return 0;
    }
    fwrite("RIFF", 1, 4, fp);
    int length = size_bytes + 44 - 8;
    fwrite(&length, 1, 4, fp);
    int is_floating = bits_per_sample < 0; if(is_floating) bits_per_sample = -bits_per_sample; //< @r-lyeh: added float support
    fwrite(is_floating ? "WAVEfmt \x10\x00\x00\x00\x03\x00":"WAVEfmt \x10\x00\x00\x00\x01\x00", 1, 14, fp); //< @r-lyeh: added float support
    fwrite(&num_channels, 1, 2, fp);
    fwrite(&sample_rate_hz, 1, 4, fp);
    int bpsec = num_channels * sample_rate_hz * bits_per_sample/8;
    fwrite(&bpsec, 1, 4, fp);
    short bpsamp = num_channels * bits_per_sample/8;
    fwrite(&bpsamp, 1, 2, fp);
    fwrite(&bits_per_sample, 1, 2, fp);
    fwrite("data", 1, 4, fp);
    fwrite(&size_bytes, 1, 4, fp);
    fwrite(samples, 1, size_bytes, fp);
    fclose(fp);
    return 1;
}

static const int jo_wav_ima_step_table[89] = { 7, 8, 9, 10, 11, 12, 13, 14, 16, 17, 19, 21, 23, 25, 28, 31, 34, 37, 41, 45, 50, 55, 60, 66, 73, 80, 88, 97, 107, 118, 130, 143, 157, 173, 190, 209, 230, 253, 279, 307, 337, 371, 408, 449, 494, 544, 598, 658, 724, 796, 876, 963, 1060, 1166, 1282, 1411, 1552, 1707, 1878, 2066, 2272, 2499, 2749, 3024, 3327, 3660, 4026, 4428, 4871, 5358, 5894, 6484, 7132, 7845, 8630, 9493, 10442, 11487, 12635, 13899, 15289, 16818, 18500, 20350, 22385, 24623, 27086, 29794, 32767 };
#define JO_WAV_INDEX_ADJUST(code) ((code & 7) < 4 ? -1 : ((code & 7) - 3) * 2)

typedef struct {
    int sample;
    int index;
    int error;         // First-order error
    int prev_error;    // Second-order error
} jo_wav_encode_channel_state_t;

#define JO_WAV_MAX_CHANNELS (8)

typedef struct {
    jo_wav_encode_channel_state_t channels[JO_WAV_MAX_CHANNELS];
    int num_channels, lookahead, quality, greedy;
    long long total_error;
} jo_wav_encode_state_t;

#define JO_WAV_CLAMP(data, min, max) ((data) < (min) ? (min) : (data) > (max) ? (max) : (data))

#ifndef JO_WAV_SAMP_CLAMP
#define JO_WAV_SAMP_CLAMP(data) ((data) < -32768 ? -32768 : (data) > 32767 ? 32767 : (data))
#endif

#define JO_WAV_ERROR_L1(data) ((data) < 0 ? -(data) : (data));
#define JO_WAV_ERROR_L2(data) ((long long)(data) * (long long)(data));

static const int jo_wav_check_order[] = {-1,+1,-8,+8,-2,+2,-4,+4,-6,+6}; // in most likely better, to least likely better order

static inline int compute_trial_delta(int code, int step) {
    // Standard IMA-ADPCM uses: diff = ((2*(code & 7) + 1) * step) >> 3
    // We add 4 (half of 8) before division to achieve rounding to nearest.
    int diff = ((2 * (code & 7) + 1) * step + 4) / 8;
    return (code & 8) ? -diff : diff;
}

static long long jo_wav_encode_adpcm_sample_r(const jo_wav_encode_state_t *state, const jo_wav_encode_channel_state_t *pchan, int num_channels, int s, const short *sample, int depth, int *out_code) {
    int delta, step, code, code2, num_check, best_code, trial_delta, mag;
    long long error, est, best, best_est;
    jo_wav_encode_channel_state_t chan;

    chan = *pchan;
    delta = s - chan.sample;
    step = jo_wav_ima_step_table[chan.index];

    // Compute the likely best answer.
    if (delta < 0) {
        mag = (-delta << 2) / step;
        code = 0x8 | (mag > 7 ? 7 : mag);
    } else {
        mag = (delta << 2) / step;
        code = mag > 7 ? 7 : mag;
    }

    // Originally computed with multiple shifts; now we use our helper:
    trial_delta = compute_trial_delta(code, step);

    chan.sample += trial_delta;
    chan.sample = JO_WAV_SAMP_CLAMP(chan.sample);
    best = JO_WAV_ERROR_L2(s - chan.sample);
    if (depth <= 0) {
        if(out_code) *out_code = code;
        return best;
    }
    
    chan.index += JO_WAV_INDEX_ADJUST(code);
    chan.index = chan.index < 0 ? 0 : chan.index > 88 ? 88 : chan.index;
#ifndef JO_WAV_NO_NOISE_SHAPING
    chan.error = chan.sample - s;
    int combined_error = (chan.error + chan.prev_error) >> 1;
    best_est = jo_wav_encode_adpcm_sample_r(state, &chan, num_channels, sample[num_channels] - combined_error, sample + num_channels, depth - 1, NULL);
#else
    best_est = jo_wav_encode_adpcm_sample_r(state, &chan, num_channels, sample[num_channels], sample + num_channels, depth - 1, NULL);
#endif
    
    // Try surrounding codes recursively by checking alternate code values.
    num_check = state->quality*2;
    num_check = JO_WAV_CLAMP(num_check, 0, (int)(sizeof(jo_wav_check_order)/sizeof(jo_wav_check_order[0])));
    best_code = code;
    for(int i = 0; i < num_check; ++i) {
        code2 = code + jo_wav_check_order[i];
        if(code2 < 0 || code2 > 15) {
            continue;
        }

        trial_delta = compute_trial_delta(code2, step);

        chan = *pchan;
        chan.sample += trial_delta;
        chan.sample = JO_WAV_SAMP_CLAMP(chan.sample);
        error = JO_WAV_ERROR_L2(s - chan.sample);
        if(!state->greedy || error < best) {
            chan.index += JO_WAV_INDEX_ADJUST(code2);
            chan.index = chan.index < 0 ? 0 : chan.index > 88 ? 88 : chan.index;
#ifndef JO_WAV_NO_NOISE_SHAPING
            chan.error = chan.sample - s;
            int combined_error_alt = (chan.error + chan.prev_error) >> 1;
            est = jo_wav_encode_adpcm_sample_r(state, &chan, num_channels, sample[num_channels] - combined_error_alt, sample + num_channels, depth - 1, NULL);
#else
            est = jo_wav_encode_adpcm_sample_r(state, &chan, num_channels, sample[num_channels], sample + num_channels, depth - 1, NULL);
#endif
            if (error + est < best + best_est) {
                best_code = code2;
                best_est = est;
                best = error;
            }
        }
    }
    
    if(out_code) *out_code = best_code;
    return best + best_est;
}

static int jo_wav_encode_adpcm_sample(jo_wav_encode_state_t *state, int ch, const short *sample, int num_samples) {
    int depth, code, step, trial_delta;
    jo_wav_encode_channel_state_t *chan = state->channels + ch;
    
    int combined_error = (chan->error + chan->prev_error) >> 1;
    int s = *sample - combined_error;

    // Add minimal dither (-1, 0, or +1) to decorrelate the quantization noise.
    s += roundf(((rand() / (float)RAND_MAX) + (rand() / (float)RAND_MAX) - 1.0f));

    depth = num_samples - 1;
    depth = depth > state->lookahead ? state->lookahead : depth;

    jo_wav_encode_adpcm_sample_r(state, chan, state->num_channels, s, sample, depth, &code);

    step = jo_wav_ima_step_table[chan->index];
    trial_delta = compute_trial_delta(code, step);

    chan->sample += trial_delta;
    chan->sample = JO_WAV_SAMP_CLAMP(chan->sample);
    chan->index += JO_WAV_INDEX_ADJUST(code);
    chan->index = chan->index < 0 ? 0 : chan->index > 88 ? 88 : chan->index;
    
#ifndef JO_WAV_NO_NOISE_SHAPING
    // Calculate current error and update the second-order shaping history.
    int current_error = chan->sample - s;
    chan->prev_error = chan->error;
    chan->error = current_error;
#endif

    state->total_error += JO_WAV_ERROR_L1(*sample - chan->sample);
    return code;
}

int jo_write_wav_adpcm_ex(const char *filename, short num_channels, int sample_rate_hz, const short *samples, int num_samples, int samples_per_block, int lookahead, int greedy, int quality) {
    FILE *fp = fopen(filename, "wb");
    if(!fp) {
        return 0;
    }

    int block_size = (samples_per_block - 1) / (num_channels ^ 3) + (num_channels * 4);
    int num_blocks = num_samples / samples_per_block;
    int total_data_bytes = num_blocks * block_size;
    int bytes_per_second = sample_rate_hz * block_size / samples_per_block;

    int leftover_samples = num_samples - num_blocks * samples_per_block;
    if (leftover_samples) {
        int last_block_samples = ((leftover_samples + 6) & ~7) + 1;
        int last_block_size = (last_block_samples - 1) / (num_channels ^ 3) + (num_channels * 4);
        total_data_bytes += last_block_size;
    }

    fwrite("RIFF", 1, 4, fp);
    int total_length = 12+20+12+8 + total_data_bytes;
    fwrite(&total_length, 1, 4, fp);
    fwrite("WAVEfmt \x14\x00\x00\x00\x11\x00", 1, 14, fp); // 0x11 == IMA_ADPCM
    fwrite(&num_channels, 1, 2, fp);
    fwrite(&sample_rate_hz, 1, 4, fp);
    fwrite(&bytes_per_second, 1, 4, fp);
    fwrite(&block_size, 1, 2, fp);
    fwrite("\x04\x00\x02\x00", 1, 4, fp); 
    fwrite(&samples_per_block, 1, 2, fp);
    fwrite("fact\x04\x00\x00\x00", 1, 8, fp);
    fwrite(&num_samples, 1, 4, fp);
    fwrite("data", 1, 4, fp);
    fwrite(&total_data_bytes, 1, 4, fp);
    
    // Now ADPCM encode the data.
    jo_wav_encode_state_t state = {0};
    state.num_channels = num_channels;
    state.lookahead = lookahead;
    state.quality = quality+1;
    state.greedy = greedy;
    for(int i = 0; i < num_samples; i += samples_per_block) {
        // initial sample
        for(int ch = 0; ch < num_channels; ++ch) {
            state.channels[ch].sample = samples[i*num_channels+ch];
            if(i == 0) for(int j = 0; j <= 88; j++) {
                if(j == 88 || state.channels[ch].sample < (jo_wav_ima_step_table[j] + jo_wav_ima_step_table[j+1] + 1) / 2) {
                    state.channels[ch].index = j;
                    break;
                }
            }
            fwrite(&state.channels[ch].sample, 2, 1, fp);
            fputc(state.channels[ch].index, fp);
            fputc(0, fp);
        }
        for(int j = i+1; j < num_samples && j < i+samples_per_block; j += 8) {
            for(int ch = 0; ch < num_channels; ++ch) {
                const short *S = samples+j*num_channels+ch;
                for(int k = 0; k < 4; ++k) {
                    unsigned char byte = jo_wav_encode_adpcm_sample(&state, ch, S, num_samples-j);
                    byte |= jo_wav_encode_adpcm_sample(&state, ch, S + num_channels, num_samples-j-1) << 4;
                    S += num_channels*2;
                    fputc(byte, fp);
                }
            }
        }
    }

#ifdef JO_WAV_TEST
    printf("Total ADPCM encoding error = %lli\n", state.total_error);
    printf("Average error per sample = %f\n", (float)state.total_error / num_samples);
#endif

    fclose(fp);
    return 1;
}

// Uses samples_per_block=2041, lookahead=0, greedy=1, quality=0
int jo_write_wav_adpcm_fastest(const char *filename, short num_channels, int sample_rate_hz, const short *samples, int num_samples) {
    return jo_write_wav_adpcm_ex(filename, num_channels, sample_rate_hz, samples, num_samples, 2041, 0, 1, 0);
}

// Uses samples_per_block=2041, lookahead=3, greedy=1, quality=4
int jo_write_wav_adpcm_fast(const char *filename, short num_channels, int sample_rate_hz, const short *samples, int num_samples) {
    return jo_write_wav_adpcm_ex(filename, num_channels, sample_rate_hz, samples, num_samples, 2041, 3, 1, 4);
}

// Uses samples_per_block=2041, lookahead=3, greedy=0, quality=4
int jo_write_wav_adpcm_slow(const char *filename, short num_channels, int sample_rate_hz, const short *samples, int num_samples) {
    return jo_write_wav_adpcm_ex(filename, num_channels, sample_rate_hz, samples, num_samples, 2041, 3, 0, 4);
}

#endif

#ifdef __cplusplus
}
#endif

#ifdef JO_WAV_TEST

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

// Test program: generates a sine sweep (logarithmic) and writes two WAV files:
// one using PCM encoding via jo_write_wav() and one using ADPCM encoding via jo_write_wav_adpcm_slow().
int main(void) {
    const int sample_rate = 44100;
    const int duration_seconds = 2;
    const int num_samples = sample_rate * duration_seconds;

    // Define the start and end frequencies for the sweep.
    double f0 = 20.0;       // Start frequency in Hz.
    double f1 = 20000.0;    // End frequency in Hz.
    
    // Precompute the sweep rate for an exponential sweep:
    double T = duration_seconds;    // Total duration
    double k = log(f1 / f0) / T;      // Exponential sweep constant
    
    // Allocate samples array.
    short *samples = (short *)malloc(num_samples * sizeof(short));
    if (!samples) {
        fprintf(stderr, "Failed to allocate memory for samples.\n");
        return 1;
    }
    
    // Generate the sine sweep.
    // The phase of an exponential sweep is computed as:
    // phase(t) = 2Ï€ * f0 * (exp(k * t) - 1) / k
    for (int i = 0; i < num_samples; ++i) {
        double t = (double)i / sample_rate;
        double phase = 2 * M_PI * f0 * (exp(k * t) - 1) / k;
        samples[i] = (short)(32767 * sin(phase));
    }

    // Write a PCM WAV file.
    printf("Writing PCM WAV file: test_pcm.wav\n");
    if (!jo_write_wav("test_pcm.wav", 1, 16, sample_rate, samples, num_samples * sizeof(short))) {
        fprintf(stderr, "Error writing PCM WAV file.\n");
        free(samples);
        return 1;
    }
    
    // Write an ADPCM WAV file using the "fast" encoder.
    printf("Writing ADPCM WAV file: test_adpcm.wav\n");
    if (!jo_write_wav_adpcm_slow("test_adpcm.wav", 1, sample_rate, samples, num_samples)) {
        fprintf(stderr, "Error writing ADPCM WAV file.\n");
        free(samples);
        return 1;
    }
    
    free(samples);
    printf("Test WAV files written successfully.\n");
    return 0;
}

#endif  // JO_WAV_TEST