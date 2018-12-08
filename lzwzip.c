#include <stdio.h>
#include <stdlib.h>
#include <string.h>

typedef struct {
  FILE* fp;
  unsigned int data;
  unsigned int bits;
  int pos;
  int num;
  unsigned char buf[4096];
} bit_buf;

bit_buf* bit_init(FILE* fp) {
  bit_buf* bib = malloc(sizeof(bit_buf));
  if (bib == NULL) return NULL;
  bib->fp = fp;
  bib->data = 0;
  bib->bits = 0;
  bib->pos = 0;
  bib->num = 0;
  return bib;
}

int bit_write(unsigned int data, unsigned int bits, bit_buf* bib) {
  bib->bits += bits;
  if (bib->bits > 32) {
    bib->bits -= 32;
    bib->data |= data >> bib->bits;
    bib->buf[bib->pos++] = bib->data >> 24;
    bib->buf[bib->pos++] = bib->data >> 16;
    bib->buf[bib->pos++] = bib->data >> 8;
    bib->buf[bib->pos++] = bib->data;
    if (bib->pos == 4096) {
      if (fwrite(bib->buf, 1, 4096, bib->fp) != 4096) return -1;
      bib->pos = 0;
    }
    bib->data = 0;
  }
  bib->data |= data << (32 - bib->bits);
  return 0;
}

int bit_flush(bit_buf* bib) {
  if (bib->bits) {
    bib->buf[bib->pos++] = bib->data >> 24;
  }
  if (bib->bits > 8) {
    bib->buf[bib->pos++] = bib->data >> 16;
  }
  if (bib->bits > 16) {
    bib->buf[bib->pos++] = bib->data >> 8;
  }
  if (bib->bits > 24) {
    bib->buf[bib->pos++] = bib->data;
  }
  if (fwrite(bib->buf, 1, bib->pos, bib->fp) != bib->pos) return -1;
  bib->bits = 0;
  bib->pos = 0;
  return 0;
}

int bit_read(unsigned int bits, bit_buf* bib) {
  int ret = 0;
  while (bits > bib->bits) {
    bits -= bib->bits;
    ret |= (bib->data & ((1 << bib->bits) - 1)) << bits;
    if (bib->pos == bib->num) {
      bib->num = fread(bib->buf, 1, 4096, bib->fp);
      if (bib->num == 0) return feof(bib->fp) ? -2 : -1;
      bib->pos = 0;
    }
    if (bib->pos + 4 <= bib->num) {
      unsigned int data = bib->buf[bib->pos];
      data <<= 8;
      data |= bib->buf[bib->pos + 1];
      data <<= 8;
      data |= bib->buf[bib->pos + 2];
      data <<= 8;
      data |= bib->buf[bib->pos + 3];
      bib->data = data;
      bib->pos += 4;
      bib->bits = 32;
    } else {
      bib->data = (unsigned int)bib->buf[bib->pos++];
      bib->bits = 8;
    }
  }
  bib->bits -= bits;
  ret |= ((bib->data >> bib->bits) & ((1 << bits) - 1));
  return ret;
}

typedef struct {
  FILE* fp;
  int pos;
  int num;
  char buf[4096];
} byte_buf;

byte_buf* byte_init(FILE* fp) {
  byte_buf* byb = malloc(sizeof(byte_buf));
  if (byb == NULL) return NULL;
  byb->fp = fp;
  byb->pos = 0;
  byb->num = 0;
  return byb;
}

int byte_read(byte_buf* byb) {
  if (byb->pos != byb->num) {
    return (unsigned char) byb->buf[byb->pos++];
  }
  byb->num = fread(byb->buf, 1, 4096, byb->fp);
  if (byb->num == 0) {
    return feof(byb->fp) ? -2 : -1;
  }
  byb->pos = 1;
  return (unsigned char) byb->buf[0];
}

int byte_write(char* seq, int size, byte_buf* byb) {
  while (size + byb->pos > 4096) {
    memcpy(byb->buf + byb->pos, seq, 4096 - byb->pos);
    if (fwrite(byb->buf, 1, 4096, byb->fp) != 4096) return -1;
    seq = &seq[4096 - byb->pos];
    size -= 4096 - byb->pos;
    byb->pos = 0;
  }
  memcpy(byb->buf + byb->pos, seq, size);
  byb->pos += size;
  return 0;
}

int byte_flush(byte_buf* byb) {
  if (byb->pos != 0) {
    if (fwrite(byb->buf, 1, byb->pos, byb->fp) != byb->pos) return -1;
    byb->pos = 0;
  }
  return 0;
}

int number_arg(char* arg) {
  size_t i;
  for (i = 0; arg[i]; ++i) {
    if (arg[i] < '0' || arg[i] > '9') return 0;
  }
  return 1;
}

void reverse(unsigned char* s, unsigned int l) {
  unsigned char* r = &s[l - 1];
  unsigned char* e = &s[l >> 1];
  while (s != e) {
    unsigned char t = *s;
    *s++ = *r;
    *r-- = t;
  }
}

int main(int argc, char* argv[]) {
  enum { OK = 0, ENC_ERROR, MEM_ERROR, IN_ERROR, OUT_ERROR, ARG_ERROR };
  bit_buf* bib = NULL;
  byte_buf* byb = NULL;
  FILE* ifp = stdin;
  FILE* ofp = stdout;
  int byte;
  int pos;
  int ppos;
  int npos;
  int arg;
  int arg_b = 0;
  int arg_d = 0;
  int arg_i = 0;
  int arg_o = 0;
  int arg_p = 0;
  int err = OK;
  unsigned long long oc = 0ull;
  unsigned long long ic = 0ull;
  unsigned long long osc = 0ull;
  unsigned long long isc = 1ull;
  unsigned int* temp;
  unsigned int* pool = NULL;
  unsigned int pool_kb = 256;
  unsigned int pool_ukb = 1;
  unsigned int bits = 9;
  unsigned int mbits = 16;
  unsigned int dir_s = 257;
  unsigned int seq_s = 1;
  unsigned char* seq = NULL;
  for (arg = 1; arg != argc; ++arg) {
    if (strcmp(argv[arg], "-h") == 0) {
      fprintf(stderr, "\
Usage: %s [OPTION] ...\n\n\
Compress or decompress files using LZW algorithm\n\n\
  -h, type this help message\n\
  -d, decompress (else execute compress)\n\
  -b, set max bits for dictionary codes (for compress)\n\
  -p, type compression ratio at standard error (for compress)\n\
  -i file, give input file (else reads from standard input)\n\
  -o file, give output file (else writes at standard output)\n\n",
              argv[0]);
      goto no_error;
    } else if (strcmp(argv[arg], "-i") == 0) {
      if (arg_i++) {
        fprintf(stderr, "%s: Already specified -i\n", argv[0]);
      } else if (++arg == argc) {
        fprintf(stderr, "%s: Missing file after -i\n", argv[0]);
      } else if ((ifp = fopen(argv[arg], "rb")) == NULL) {
        fprintf(stderr, "%s: %s: Cannot open input file\n", argv[0], argv[arg]);
      } else {
        continue;
      }
      goto arg_error;
    } else if (strcmp(argv[arg], "-o") == 0) {
      if (arg_o++) {
        fprintf(stderr, "%s: Already specified -o\n", argv[0]);
      } else if (++arg == argc) {
        fprintf(stderr, "%s: Missing file after -o\n", argv[0]);
      } else if ((ofp = fopen(argv[arg], "wb")) == NULL) {
        fprintf(stderr, "%s: %s: Cannot open output file\n", argv[0],
                argv[arg]);
      } else {
        continue;
      }
      goto arg_error;
    } else if (strcmp(argv[arg], "-b") == 0) {
      if (arg_b++) {
        fprintf(stderr, "%s: Already specified -b\n", argv[0]);
      } else if (arg_d) {
        fprintf(stderr, "%s: Cannot use -b with -d\n", argv[0]);
      } else if (++arg == argc) {
        fprintf(stderr, "%s: Missing number after -b\n", argv[0]);
      } else if (!number_arg(argv[arg])) {
        fprintf(stderr, "%s: Not a number after -b\n", argv[0]);
      } else if ((mbits = (unsigned)atoi(argv[arg])) < 9 || mbits > 16) {
        fprintf(stderr, "%s: %d: Invalid max bits\n", argv[0], mbits);
      } else {
        continue;
      }
      goto arg_error;
    } else if (strcmp(argv[arg], "-d") == 0) {
      if (arg_d++) {
        fprintf(stderr, "%s: Already specified -d\n", argv[0]);
      } else if (arg_b) {
        fprintf(stderr, "%s: Cannot use -d with -b\n", argv[0]);
      } else if (arg_p) {
        fprintf(stderr, "%s: Cannot use -d with -p\n", argv[0]);
      } else {
        continue;
      }
      goto arg_error;
    } else if (strcmp(argv[arg], "-p") == 0) {
      if (arg_p++) {
        fprintf(stderr, "%s: Already specified -p\n", argv[0]);
      } else if (arg_d) {
        fprintf(stderr, "%s: Cannot use -p with -d\n", argv[0]);
      } else {
        continue;
      }
      goto arg_error;
    } else {
      fprintf(stderr, "%s: Invalid arguments\n", argv[0]);
      fprintf(stderr, "For help, type %s -h\n", argv[0]);
      goto arg_error;
    }
  }
  if (arg_d == 0) {
    byb = byte_init(ifp);
    if (byb == NULL) goto mem_error;
    byte = byte_read(byb);
    if (byte == -1) goto no_error;
    if (byte < -1) goto in_error;
    bib = bit_init(ofp);
    pool = malloc(pool_kb << 10);
    if (bib == NULL || pool == NULL) goto mem_error;
    if (bit_write(mbits - 9, 3, bib) == -1) goto out_error;
    for (pos = 0; pos != 256; ++pos) {
      pool[pos] = pos;
    }
    pos = byte;
    while ((byte = byte_read(byb)) >= 0) {
      ++isc;
      npos = (pool[pos] >> 16) << 8;
      if (npos && pool[npos + byte]) {
        pos = npos + byte;
        continue;
      }
      if (bit_write((pool[pos] & 0xFFFF), bits, bib) == -1) goto out_error;
      osc += bits;
      if (dir_s == (1 << bits)) {
        bits += (bits != mbits);
        if ((osc >> 3) > isc) {
          for (pos = 0; pos != 256; ++pos) {
            pool[pos] = pos;
          }
          if (bit_write(256, bits, bib)) goto out_error;
          oc += osc + bits;
          ic += isc;
          isc = osc = 0;
          pool_ukb = 1;
          dir_s = 257;
          bits = 9;
        }
      }
      if (dir_s != (1 << mbits) && isc) {
        if (npos != 0) {
          pool[npos + byte] = dir_s++;
        } else {
          if (pool_ukb == pool_kb) {
            pool_kb <<= 1;
            temp = realloc(pool, pool_kb << 10);
            if (temp == NULL) goto mem_error;
            pool = temp;
          }
          memset(&pool[pool_ukb << 8], 0, 1024);
          pool[pos] |= (pool_ukb << 16);
          pool[((pool_ukb++) << 8) + byte] = dir_s++;
        }
      }
      pos = byte;
    }
    if (byte == -1) goto in_error;
    if (bit_write(pool[pos] & 0xFFFF, bits, bib) == -1) goto out_error;
    oc = (oc + osc + bits + 10) >> 3;
    ic += isc;
    if (bit_flush(bib) == -1) goto out_error;
  } else {
    if ((bib = bit_init(ifp)) == NULL) goto mem_error;
    mbits = bit_read(3, bib);
    if (mbits == -1) goto no_error;
    if (mbits < -1) goto in_error;
    mbits += 9;
    byb = byte_init(ofp);
    seq = malloc(1 << mbits);
    pool = malloc(1 << (mbits + 2));
    if (byb == NULL || seq == NULL || pool == NULL) goto mem_error;
    for (pos = 0; pos != 256; ++pos) {
      pool[pos] = pos;
    }
    do {
      ppos = bit_read(bits, bib);
    } while (ppos == 256);
    if (ppos == -1) goto enc_error;
    seq[0] = ppos;
    if (byte_write((char*) seq, 1, byb) == -1) goto out_error;
    while ((npos = bit_read(bits, bib)) >= 0) {
      if (npos > (int) dir_s) goto enc_error;
      if (npos == 256) {
        bits = 9;
        do {
          ppos = bit_read(bits, bib);
        } while (ppos == 256);
        if (ppos == -1) break;
        seq[0] = ppos;
        if (byte_write((char*) seq, 1, byb) == -1) goto out_error;
        dir_s = 257;
        seq_s = 1;
        continue;
      }
      if (npos == (int) dir_s) {
        seq[seq_s++] = seq[0];
      } else {
        pos = npos;
        seq[0] = pool[npos];
        seq_s = 1;
        while (pos > 256) {
          pos = pool[pos] >> 8;
          seq[seq_s++] = pool[pos];
        }
        reverse(seq, seq_s);
      }
      if (byte_write((char*) seq, seq_s, byb) == -1) goto out_error;
      if (dir_s != (1 << mbits)) {
        pool[dir_s++] = (ppos << 8) | seq[0];
        if (bits != mbits && dir_s == (1 << bits)) ++bits;
      }
      ppos = npos;
    }
    if (npos == -1) goto in_error;
    if (byte_flush(byb) == -1) goto out_error;
  }
no_error:
  if (fclose(ifp) == EOF && err == OK) err = IN_ERROR;
  if (fflush(ofp) == EOF && err == OK) err = OUT_ERROR;
  if (fclose(ofp) == EOF && err == OK) err = OUT_ERROR;
  free(seq);
  free(pool);
  free(byb);
  free(bib);
  if (err == ENC_ERROR) {
    fprintf(stderr, "%s: Invalid input file encoding\n", argv[0]);
  } else if (err == MEM_ERROR) {
    fprintf(stderr, "%s: Not enough memory\n", argv[0]);
  } else if (err == IN_ERROR) {
    fprintf(stderr, "%s: Error while reading input file\n", argv[0]);
  } else if (err == OUT_ERROR) {
    fprintf(stderr, "%s: Error while writing output file\n", argv[0]);
  } else if (err == OK && arg_p) {
    double ratio = 100 - 100.0 * oc / ic;
    if (oc <= ic) {
      fprintf(stderr, "Achieved compression of %.2f%%\n", ratio);
    } else {
      fprintf(stderr, "Used additional %.2f%%\n", -ratio);
    }
  }
  return err;
enc_error:
  err = ENC_ERROR;
  goto no_error;
mem_error:
  err = MEM_ERROR;
  goto no_error;
in_error:
  err = IN_ERROR;
  goto no_error;
out_error:
  err = OUT_ERROR;
  goto no_error;
arg_error:
  err = ARG_ERROR;
  goto no_error;
}
