#include <stdlib.h>
#include <png.h>
#include <zlib.h>
#include <assert.h>

#include "bitmap.h"

void write_png(FILE* file, int width, int height, unsigned char* buffer) {
	png_bytep* row_pointers=(png_bytep*)malloc(sizeof(png_bytep)*height); assert(row_pointers);
	for(int c=0;c<height;c++) row_pointers[c]=buffer+4*width*(height-c-1);
    /* initialize stuff */
    png_structp png_ptr=png_create_write_struct(PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);
    assert(png_ptr);

    png_infop info_ptr = png_create_info_struct(png_ptr);
	assert(info_ptr);

    assert(!setjmp(png_jmpbuf(png_ptr)));

    png_init_io(png_ptr, file);

    /* write header */
    assert(!setjmp(png_jmpbuf(png_ptr)));

    png_set_IHDR(png_ptr, info_ptr, width, height,
        8, PNG_COLOR_TYPE_RGB_ALPHA, PNG_INTERLACE_NONE,
        PNG_COMPRESSION_TYPE_BASE, PNG_FILTER_TYPE_BASE);

    png_write_info(png_ptr, info_ptr);

    /* write bytes */
    assert(!setjmp(png_jmpbuf(png_ptr)));

    png_write_image(png_ptr, row_pointers);
    /* end write */
    assert(!setjmp(png_jmpbuf(png_ptr)));

    png_write_end(png_ptr, NULL);
    free(row_pointers);
}