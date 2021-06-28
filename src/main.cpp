#include <stdio.h>
#include <setjmp.h>
#include <boost/scoped_ptr.hpp>
#include <boost/numeric/conversion/cast.hpp>
#include <vector>
#include <cstdlib>
#include <iostream>

#include <fstream>

#include <string.h> // for memcpy

extern "C" {
    #include <jpeglib.h>
    #include "jerror.h"
}

/* File is compiled inside one of the *bit directories*/
// #include "../jpeg_utils.h"

#include <assert.h>

#define MEX_JPEG_MODE    int
#define MEX_JPEG_LOSSY     0
#define MEX_JPEG_LOSSLESS  1

#define UCHAR UINT8_T

#define MEX_JPEG_MODE    int
#define MEX_JPEG_LOSSY     0
#define MEX_JPEG_LOSSLESS  1

#define UCHAR UINT8_T

// --------------------------------------------------------------------------------------------------------------------------------------------
#define SIZEOF(object)	((size_t) sizeof(object))
#define JFWRITE(file,buf,sizeofbuf)  \
  ((size_t) fwrite((const void *) (buf), (size_t) 1, (size_t) (sizeofbuf), (file)))
#define MEMCOPY(dest,src,size)	memcpy((void *)(dest), (const void *)(src), (size_t)(size))

#ifndef HAVE_STDLIB_H		/* <stdlib.h> should declare malloc(),free() */
extern void * malloc JPP((size_t size));
extern void free JPP((void *ptr));
#endif


/* Expanded data destination object for stdio output */

typedef struct {
  struct jpeg_destination_mgr pub; /* public fields */

  FILE * outfile;		/* target stream */
  JOCTET * buffer;		/* start of buffer */
} my_destination_mgr;

typedef my_destination_mgr * my_dest_ptr;

#define OUTPUT_BUF_SIZE  4096	/* choose an efficiently fwrite'able size */


/* Expanded data destination object for memory output */

typedef struct {
  struct jpeg_destination_mgr pub; /* public fields */

  unsigned char ** outbuffer;	/* target buffer */
  unsigned long * outsize;
  unsigned char * newbuffer;	/* newly allocated buffer */
  JOCTET * buffer;		/* start of buffer */
  size_t bufsize;
} my_mem_destination_mgr;

typedef my_mem_destination_mgr * my_mem_dest_ptr;


/*
 * Initialize destination --- called by jpeg_start_compress
 * before any data is actually written.
 */

METHODDEF(void)
init_destination (j_compress_ptr cinfo)
{
  my_dest_ptr dest = (my_dest_ptr) cinfo->dest;

  /* Allocate the output buffer --- it will be released when done with image */
  dest->buffer = (JOCTET *)
      (*cinfo->mem->alloc_small) ((j_common_ptr) cinfo, JPOOL_IMAGE,
				  OUTPUT_BUF_SIZE * SIZEOF(JOCTET));

  dest->pub.next_output_byte = dest->buffer;
  dest->pub.free_in_buffer = OUTPUT_BUF_SIZE;
}

METHODDEF(void)
init_mem_destination (j_compress_ptr cinfo)
{
  /* no work necessary here */
}


/*
 * Empty the output buffer --- called whenever buffer fills up.
 *
 * In typical applications, this should write the entire output buffer
 * (ignoring the current state of next_output_byte & free_in_buffer),
 * reset the pointer & count to the start of the buffer, and return TRUE
 * indicating that the buffer has been dumped.
 *
 * In applications that need to be able to suspend compression due to output
 * overrun, a FALSE return indicates that the buffer cannot be emptied now.
 * In this situation, the compressor will return to its caller (possibly with
 * an indication that it has not accepted all the supplied scanlines).  The
 * application should resume compression after it has made more room in the
 * output buffer.  Note that there are substantial restrictions on the use of
 * suspension --- see the documentation.
 *
 * When suspending, the compressor will back up to a convenient restart point
 * (typically the start of the current MCU). next_output_byte & free_in_buffer
 * indicate where the restart point will be if the current call returns FALSE.
 * Data beyond this point will be regenerated after resumption, so do not
 * write it out when emptying the buffer externally.
 */

METHODDEF(boolean)
empty_output_buffer (j_compress_ptr cinfo)
{
  my_dest_ptr dest = (my_dest_ptr) cinfo->dest;

  if (JFWRITE(dest->outfile, dest->buffer, OUTPUT_BUF_SIZE) !=
      (size_t) OUTPUT_BUF_SIZE)
    ERREXIT(cinfo, JERR_FILE_WRITE);

  dest->pub.next_output_byte = dest->buffer;
  dest->pub.free_in_buffer = OUTPUT_BUF_SIZE;

  return TRUE;
}

METHODDEF(boolean)
empty_mem_output_buffer (j_compress_ptr cinfo)
{
  size_t nextsize;
  JOCTET * nextbuffer;
  my_mem_dest_ptr dest = (my_mem_dest_ptr) cinfo->dest;

  /* Try to allocate new buffer with double size */
  nextsize = dest->bufsize * 2;
  nextbuffer = (JOCTET*)(malloc(nextsize));

  if (nextbuffer == NULL)
    ERREXIT1(cinfo, JERR_OUT_OF_MEMORY, 10);

  MEMCOPY(nextbuffer, dest->buffer, dest->bufsize);

  if (dest->newbuffer != NULL)
    free(dest->newbuffer);

  dest->newbuffer = nextbuffer;

  dest->pub.next_output_byte = nextbuffer + dest->bufsize;
  dest->pub.free_in_buffer = dest->bufsize;

  dest->buffer = nextbuffer;
  dest->bufsize = nextsize;

  return TRUE;
}


/*
 * Terminate destination --- called by jpeg_finish_compress
 * after all data has been written.  Usually needs to flush buffer.
 *
 * NB: *not* called by jpeg_abort or jpeg_destroy; surrounding
 * application must deal with any cleanup that should happen even
 * for error exit.
 */

METHODDEF(void)
term_destination (j_compress_ptr cinfo)
{
  my_dest_ptr dest = (my_dest_ptr) cinfo->dest;
  size_t datacount = OUTPUT_BUF_SIZE - dest->pub.free_in_buffer;

  /* Write any data remaining in the buffer */
  if (datacount > 0) {
    if (JFWRITE(dest->outfile, dest->buffer, datacount) != datacount)
      ERREXIT(cinfo, JERR_FILE_WRITE);
  }
  fflush(dest->outfile);
  /* Make sure we wrote the output file OK */
  if (ferror(dest->outfile))
    ERREXIT(cinfo, JERR_FILE_WRITE);
}

METHODDEF(void)
term_mem_destination (j_compress_ptr cinfo)
{
  my_mem_dest_ptr dest = (my_mem_dest_ptr) cinfo->dest;

  *dest->outbuffer = dest->buffer;
  *dest->outsize = dest->bufsize - dest->pub.free_in_buffer;
}


/*
 * Prepare for output to a stdio stream.
 * The caller must have already opened the stream, and is responsible
 * for closing it after finishing compression.
 */

GLOBAL(void)
jpeg_stdio_dest (j_compress_ptr cinfo, FILE * outfile)
{
  my_dest_ptr dest;

  /* The destination object is made permanent so that multiple JPEG images
   * can be written to the same file without re-executing jpeg_stdio_dest.
   * This makes it dangerous to use this manager and a different destination
   * manager serially with the same JPEG object, because their private object
   * sizes may be different.  Caveat programmer.
   */
  if (cinfo->dest == NULL) {	/* first time for this JPEG object? */
    cinfo->dest = (struct jpeg_destination_mgr *)
      (*cinfo->mem->alloc_small) ((j_common_ptr) cinfo, JPOOL_PERMANENT,
				  SIZEOF(my_destination_mgr));
  }

  dest = (my_dest_ptr) cinfo->dest;
  dest->pub.init_destination = init_destination;
  dest->pub.empty_output_buffer = empty_output_buffer;
  dest->pub.term_destination = term_destination;
  dest->outfile = outfile;
}


/*
 * Prepare for output to a memory buffer.
 * The caller may supply an own initial buffer with appropriate size.
 * Otherwise, or when the actual data output exceeds the given size,
 * the library adapts the buffer size as necessary.
 * The standard library functions malloc/free are used for allocating
 * larger memory, so the buffer is available to the application after
 * finishing compression, and then the application is responsible for
 * freeing the requested memory.
 */

GLOBAL(void)
jpeg_mem_dest (j_compress_ptr cinfo,
	       unsigned char ** outbuffer, unsigned long * outsize)
{
  my_mem_dest_ptr dest;

  if (outbuffer == NULL || outsize == NULL)	/* sanity check */
    ERREXIT(cinfo, JERR_BUFFER_SIZE);

  /* The destination object is made permanent so that multiple JPEG images
   * can be written to the same buffer without re-executing jpeg_mem_dest.
   */
  if (cinfo->dest == NULL) {	/* first time for this JPEG object? */
    cinfo->dest = (struct jpeg_destination_mgr *)
      (*cinfo->mem->alloc_small) ((j_common_ptr) cinfo, JPOOL_PERMANENT,
				  SIZEOF(my_mem_destination_mgr));
  }

  dest = (my_mem_dest_ptr) cinfo->dest;
  dest->pub.init_destination = init_mem_destination;
  dest->pub.empty_output_buffer = empty_mem_output_buffer;
  dest->pub.term_destination = term_mem_destination;
  dest->outbuffer = outbuffer;
  dest->outsize = outsize;
  dest->newbuffer = NULL;

  if (*outbuffer == NULL || *outsize == 0) {
    /* Allocate initial buffer */
    dest->newbuffer = *outbuffer = (JOCTET*)(malloc(OUTPUT_BUF_SIZE));
    if (dest->newbuffer == NULL)
      ERREXIT1(cinfo, JERR_OUT_OF_MEMORY, 10);
    *outsize = OUTPUT_BUF_SIZE;
  }

  dest->pub.next_output_byte = dest->buffer = *outbuffer;
  dest->pub.free_in_buffer = dest->bufsize = *outsize;
}

// --------------------------------------------------------------------------------------------------------------------------------------------

namespace jpegutils {
    static char err_warn_buffer[JMSG_LENGTH_MAX];

    /* 
     * Used as the error manager by the jpeg library in read_jpeg.cpp and write_jpeg.cpp
     */
    struct my_error_mgr {
      struct jpeg_error_mgr pub;    /* "public" fields */
      jmp_buf setjmp_buffer;    /* for return to caller */
      int     read_done;        /* If reading is done then only warn not error */
    };
    typedef struct my_error_mgr* my_error_ptr;
}

static void
my_output_message (j_common_ptr cinfo)
{

  /* Create the message */
  (*cinfo->err->format_message) (cinfo, jpegutils::err_warn_buffer);

  /* cinfo->err really points to a my_error_mgr struct, so coerce pointer */
  jpegutils::my_error_ptr myerr = reinterpret_cast<jpegutils::my_error_ptr>(cinfo->err);

  /* We cannot continue execution of jpeg_finish_decompress during cleanup after
   * we reach one of these errors. The behavior of the JPEG libaray is not
   * defined if we continue execution and it may crash MATLAB.
   *
   * This exception is caught in the mexFunction() of read_jpeg.cpp
   * (g987471), (g135392)
  */
  if( cinfo->is_decompressor )
  {
      if (myerr->read_done)
      {
          if (cinfo->err->msg_code >= JERR_ARITH_NOTIMPL &&
              cinfo->err->msg_code <= JERR_XMS_WRITE)
          {
              throw std::runtime_error("Fatal jpeg_finish_decompress error");
          }
          else
          {
              return;
          }
      }
  }
  else
  {
      return;
  }
  
  /*
   * If any trace or warning messages are received, 
   * continue on treating the condition as a warning.
   * */

  /* Send it to stderr, adding a newline */
  // mxWarningMsgId(MATLAB::imagesci::jpg::libraryMessage(BITS_IN_JSAMPLE, jpegutils::err_warn_buffer));

}

/* 
 *  These replace specific routines in the jpeg library's
 *  jerror.c module. 
 */

enum mxClassID {
    mxUINT8_CLASS = 8,
    mxUINT16_CLASS = 16
};

struct Meta_Data {
    std::size_t ndims;
    J_COLOR_SPACE imageType;
    const size_t *size;
    mxClassID  inClass;
    volatile int quality;
    volatile MEX_JPEG_MODE compMode;
    const char *commentArray;
};

static void my_error_exit (j_common_ptr cinfo);

// static MEX_JPEG_MODE getCompMode(std::vector<unsigned char> props);

/*static J_COLOR_SPACE getImageType(size_t ndims, 
                               const size_t *size, 
                               mxClassID inClass);*/

/* static void WriteComment(j_compress_ptr cinfoPtr, 
                         std::vector<char> commentArray); */

/*static void WriteComment(j_compress_ptr cinfoPtr, 
                         const char *commentArray);*/

static void WriteGRAYFromUint8(j_compress_ptr cinfoPtr, 
                               uint8_t *uint8Data, 
                               const size_t *size);

static void WriteImageData(j_compress_ptr cinfoPtr, 
                           uint8_t *inputArray,
                           J_COLOR_SPACE imageType,
                           const size_t *size);

static void WriteImageData(j_compress_ptr cinfoPtr, 
                           uint16_t *inputArray,
                           J_COLOR_SPACE imageType,
                           const size_t *size);

static void WriteRGBFromUint8(j_compress_ptr cinfoPtr, 
                              uint8_t *uint8Data, 
                              const size_t *size);

static void WriteGRAYFromUint16(j_compress_ptr cinfoPtr, 
                               uint16_t *uint16Data, 
                               const size_t *size);

static void WriteRGBFromUint16(j_compress_ptr cinfoPtr, 
                              uint16_t *uint16Data, 
                              const size_t *size);

class JpegWriter {
    private:

    FILE *outfile;

    Meta_Data writer_meta_data;

    struct jpeg_compress_struct cinfo;

    struct jpegutils::my_error_mgr jerr;

    jpeg_destination_mgr _jdst;

    void file_init() {
        cinfo.err = jpeg_std_error(&jerr.pub);
        jerr.pub.output_message = my_output_message;
        jerr.pub.error_exit = my_error_exit;
        if(setjmp(jerr.setjmp_buffer))
        {
            /* If we get here, the JPEG code has signaled an error.
            * We need to clean up the JPEG object, close the input file,
            * and return.
            */
            jpeg_destroy_compress(&cinfo);
            if( outfile != nullptr )
            {
                fclose(outfile);
            }

            // mxErrMsgId(MATLAB::imagesci::jpg::libraryMessage(BITS_IN_JSAMPLE, jpegutils::err_warn_buffer));
            
            return;
        }

        jpeg_create_compress(&cinfo);

        cinfo.image_width  = (int) writer_meta_data.size[1];
        cinfo.image_height = (int) writer_meta_data.size[0]; 

        cinfo.in_color_space = writer_meta_data.imageType;
        if ( writer_meta_data.imageType == JCS_RGB ) {
            cinfo.input_components = 3;
        } else {
            // JCS_GRAYSCALE
            cinfo.input_components = 1;
        }


        /* Set default compression parameters. */
        jpeg_set_defaults(&cinfo);
        jpeg_set_quality(&cinfo, writer_meta_data.quality, TRUE);

        
        // outfile = fopen(filename.c_str(), "wb");
        // mxAssert((outfile != NULL), "Could not open file for writing.");

        jpeg_stdio_dest(&cinfo, outfile);

        if (writer_meta_data.compMode == MEX_JPEG_LOSSLESS) {
            jpeg_simple_lossless(&cinfo, 1, 0);
        }

        jpeg_start_compress(&cinfo, TRUE);
    }

    void buffer_init() {
        cinfo.err = jpeg_std_error(&jerr.pub);
        jerr.pub.output_message = my_output_message;
        jerr.pub.error_exit = my_error_exit;
        if(setjmp(jerr.setjmp_buffer))
        {
            /* If we get here, the JPEG code has signaled an error.
            * We need to clean up the JPEG object, close the input file,
            * and return.
            */
            jpeg_destroy_compress(&cinfo);
            if( outfile != nullptr )
            {
                fclose(outfile);
            }

            // mxErrMsgId(MATLAB::imagesci::jpg::libraryMessage(BITS_IN_JSAMPLE, jpegutils::err_warn_buffer));
            
            return;
        }

        jpeg_create_compress(&cinfo);

        cinfo.image_width  = (int) writer_meta_data.size[1];
        cinfo.image_height = (int) writer_meta_data.size[0]; 

        cinfo.in_color_space = writer_meta_data.imageType;
        if ( writer_meta_data.imageType == JCS_RGB ) {
            cinfo.input_components = 3;
        } else {
            // JCS_GRAYSCALE
            cinfo.input_components = 1;
        }


        /* Set default compression parameters. */
        jpeg_set_defaults(&cinfo);
        jpeg_set_quality(&cinfo, writer_meta_data.quality, TRUE);

        

        if (writer_meta_data.compMode == MEX_JPEG_LOSSLESS) {
            jpeg_simple_lossless(&cinfo, 1, 0);
        }

        jpeg_start_compress(&cinfo, TRUE);
    }

    void clean_up() {
        jpeg_finish_compress(&cinfo); fclose(outfile);
        jpeg_destroy_compress(&cinfo);
    }

    public:

    JpegWriter(std::string filename, Meta_Data writer_meta_data): 
    outfile(fopen(filename.c_str(), "wb")),
    writer_meta_data{
        writer_meta_data.ndims,
        writer_meta_data.imageType,
        writer_meta_data.size,
        writer_meta_data.inClass,
        writer_meta_data.quality,
        writer_meta_data.compMode,
        writer_meta_data.commentArray
    } {
        file_init();
    }

    JpegWriter(Meta_Data writer_meta_data): 
    outfile(nullptr),
    writer_meta_data{
        writer_meta_data.ndims,
        writer_meta_data.imageType,
        writer_meta_data.size,
        writer_meta_data.inClass,
        writer_meta_data.quality,
        writer_meta_data.compMode,
        writer_meta_data.commentArray
    } {
        buffer_init();
    }

    void write(uint8_t *inputArray) {
        // WriteComment(&cinfo, writer_meta_data.commentArray);
        WriteImageData(&cinfo, inputArray, writer_meta_data.imageType, writer_meta_data.size);
    }

    void write(uint16_t *inputArray) {
        // WriteComment(&cinfo, writer_meta_data.commentArray);
        WriteImageData(&cinfo, inputArray, writer_meta_data.imageType, writer_meta_data.size);
    }

    ~JpegWriter() {
        clean_up();
}

};


static void 
WriteRGBFromUint8(j_compress_ptr cinfoPtr, 
                  uint8_t *uint8Data, 
                  const size_t * size)
{
    JSAMPARRAY buffer;
    
    size_t height = size[0];
    size_t width  = size[1];
    
    JDIMENSION row_stride = boost::numeric_cast<JDIMENSION>(width * 3);    
    buffer = (*cinfoPtr->mem->alloc_sarray)
        ((j_common_ptr) cinfoPtr, JPOOL_IMAGE, row_stride, 1);
    
    /* jpeg_write_scanlines expects an array of pointers to scanlines.
     * Here the array is only one element long, but you could pass
     * more than one scanline at a time if that's more convenient.
     */
    
    while (cinfoPtr->next_scanline < cinfoPtr->image_height) {
        /* construct a buffer which contains the next scanline with data in
         * RGBRGBRGB... format */
        JDIMENSION i = cinfoPtr->next_scanline;
        for(JDIMENSION j = 0; j < cinfoPtr->image_width; j++)
        {
            buffer[0][3*j]   = uint8Data[i + (j*height)]; /* Red */
            buffer[0][3*j+1] = uint8Data[i + (j*height) + (height * width)]; /* Green */
            buffer[0][3*j+2] = uint8Data[i + (j*height) + (2 * height * width)]; /* Blue */
        }
        (void) jpeg_write_scanlines(cinfoPtr, buffer, 1);
    } 
}  
  


static void 
WriteGRAYFromUint8(j_compress_ptr cinfoPtr, 
                   uint8_t *uint8Data, 
                   const size_t *size)
{
    JSAMPARRAY buffer;
    
    size_t height = size[0];
    size_t width  = size[1];
    
    JDIMENSION row_stride = boost::numeric_cast<JDIMENSION>(width);
    buffer = (*cinfoPtr->mem->alloc_sarray)
        ((j_common_ptr) cinfoPtr, JPOOL_IMAGE, row_stride, 1);
    
    /* jpeg_write_scanlines expects an array of pointers to scanlines.
     * Here the array is only one element long, but you could pass
     * more than one scanline at a time if that's more convenient.
     */
    
    while (cinfoPtr->next_scanline < cinfoPtr->image_height) {
        /* construct a buffer which contains the next scanline with data in
         * RGBRGBRGB... format */
        JDIMENSION i = cinfoPtr->next_scanline;
        for(JDIMENSION j = 0; j < cinfoPtr->image_width; j++)
        {
            buffer[0][j]   = uint8Data[i + (j*height)]; 
        }
        (void) jpeg_write_scanlines(cinfoPtr, buffer, 1);
    } 
}  



static void 
WriteRGBFromUint16(j_compress_ptr cinfoPtr, 
                   uint16_t *uint16Data, 
                   const size_t *size)
{
    JSAMPARRAY buffer;
    
    size_t height = size[0];
    size_t width = size[1];
    
    JDIMENSION row_stride = boost::numeric_cast<JDIMENSION>(width * 3);
    buffer = (*cinfoPtr->mem->alloc_sarray)
        ((j_common_ptr) cinfoPtr, JPOOL_IMAGE, row_stride, 1);
    
    /* jpeg_write_scanlines expects an array of pointers to scanlines.
     * Here the array is only one element long, but you could pass
     * more than one scanline at a time if that's more convenient.
     */
    
    while (cinfoPtr->next_scanline < cinfoPtr->image_height) {
        /* construct a buffer which contains the next scanline with data in
         * RGBRGBRGB... format */
        JDIMENSION i = cinfoPtr->next_scanline;
        for(JDIMENSION j = 0; j < cinfoPtr->image_width; j++)
        {
            buffer[0][3*j]   = (JSAMPLE) uint16Data[i + (j*height)]; /* Red */
            buffer[0][3*j+1] = (JSAMPLE) uint16Data[i + (j*height) + (height * width)]; /* Green */
            buffer[0][3*j+2] = (JSAMPLE) uint16Data[i + (j*height) + (2 * height * width)]; /* Blue */
        }
        (void) jpeg_write_scanlines(cinfoPtr, buffer, 1);
    } 
}  



static void 
WriteGRAYFromUint16(j_compress_ptr cinfoPtr, 
                  uint16_t *uint16Data, 
                  const size_t *size)
{
    JSAMPARRAY buffer;
    
    size_t height = size[0];
    size_t width  = size[1];
    
    JDIMENSION row_stride = boost::numeric_cast<JDIMENSION>(width);
    buffer = (*cinfoPtr->mem->alloc_sarray)
        ((j_common_ptr) cinfoPtr, JPOOL_IMAGE, row_stride, 1);
    
    /* jpeg_write_scanlines expects an array of pointers to scanlines.
     * Here the array is only one element long, but you could pass
     * more than one scanline at a time if that's more convenient.
     */
    
    while (cinfoPtr->next_scanline < cinfoPtr->image_height) {
        /* construct a buffer which contains the next scanline with data in
         * RGBRGBRGB... format */
        JDIMENSION i = cinfoPtr->next_scanline;
        for(JDIMENSION j = 0; j < cinfoPtr->image_width; j++)
        {
            buffer[0][j]   = (JSAMPLE) uint16Data[i + (j*height)]; 
        }
        (void) jpeg_write_scanlines(cinfoPtr, buffer, 1);
    } 
}  



/*
 * Here's the routine that will replace the standard error_exit method:
 */

static void
my_error_exit (j_common_ptr cinfo)
{
    /* cinfo->err really points to a my_error_mgr struct, so coerce pointer */
    jpegutils::my_error_ptr myerr = (jpegutils::my_error_ptr) cinfo->err;
    
    /* Always display the message. */
    /* We could postpone this until after returning, if we chose. */
    (*cinfo->err->output_message) (cinfo);
    
    /* Return control to the setjmp point */
    longjmp(myerr->setjmp_buffer, 1);
}

static void WriteImageData(j_compress_ptr cinfoPtr, 
                           uint8_t *inputArray,
                           J_COLOR_SPACE imageType,
                           const size_t *size)
{



    if ( imageType == JCS_RGB ) {
        WriteRGBFromUint8(cinfoPtr, (uint8_t *) inputArray, size);
    } else {
        // JCS_GRAYSCALE:
        WriteGRAYFromUint8(cinfoPtr, (uint8_t *) inputArray, size);
    }

}

static void WriteImageData(j_compress_ptr cinfoPtr, 
                           uint16_t *inputArray,
                           J_COLOR_SPACE imageType,
                           const size_t *size) {

    if ( imageType == JCS_RGB ) {
        WriteRGBFromUint16(cinfoPtr, (uint16_t *) inputArray, size);
    } else {
        // JCS_GRAYSCALE:
        WriteGRAYFromUint16(cinfoPtr, (uint16_t *) inputArray, size);
    }

}


// static J_COLOR_SPACE getImageType(size_t ndims, const size_t *size, mxClassID inClass)
// {
//     J_COLOR_SPACE imgtype(JCS_RGB);

//     if( (ndims == 3 && size[2] == 3) && inClass == mxUINT8_CLASS ) {

//         imgtype = JCS_RGB;

//     } else if( (ndims == 3 && size[2] == 3) && inClass == mxUINT16_CLASS ) {

//         /*
//          * Should hold for both 12bit and 16bit
//          * */
//         imgtype = JCS_RGB;

//     } else if(ndims == 2 && inClass == mxUINT8_CLASS ) {

//         imgtype = JCS_GRAYSCALE;

//     } else if(ndims == 2 && inClass == mxUINT16_CLASS ) {

//         /*
//          * Should hold for both 12bit and 16bit
//          * */
//         imgtype = JCS_GRAYSCALE;

//     } else {

//         // mxErrMsgId(MATLAB::imagesci::jpg::invalidImage());

//     }

//     return(imgtype);

// }

int main() {
    Meta_Data writer_meta_data;

    writer_meta_data.ndims = 3;
    writer_meta_data.imageType = JCS_RGB;
    size_t *size = new size_t(3 * sizeof(size_t));
        size[0] = 87;
        size[1] = 133;
        size[2] = 3;
    writer_meta_data.size = size;
    writer_meta_data.inClass = mxUINT8_CLASS;
    writer_meta_data.quality = 75;
    writer_meta_data.compMode = MEX_JPEG_LOSSY;
    writer_meta_data.commentArray = "";

    std::string infilename = "Resources/rosew1.data";
    
    FILE * pFile = std::fopen(infilename.c_str(), "rb");
    if (pFile == NULL) {
        fputs("File error", stderr);
        exit(1);
    }

    // obtain file size:
    fseek(pFile, 0, SEEK_END);
    size_t lSize = ftell(pFile);
    rewind(pFile);

    // allocate memory to contain the whole file:
    std::vector<uint8_t> p_buffer(size[0] * size[1] * size[2] * sizeof(uint8_t));
    if (p_buffer.empty()) {
        fputs("Memory error", stderr);
        exit(2);
    }

    // copy the file contents into the buffer:
    size_t result = fread(&p_buffer[0], sizeof(uint8_t), lSize, pFile);
    if (result != lSize) {
        fputs("Reading error", stderr);
        exit(3);
    }

    fclose(pFile);

    uint8_t* uint8Data = &p_buffer[0];

    // std::string outfilename = "rosew1_copy1.jpg";
    // JpegWriter writer1(outfilename, writer_meta_data);
    // writer1.write(&p_buffer[0]);

    // // std::vector<unsigned char> outbuffer(400000);
    // // int outsize = 0;
    // JpegWriter writer2(writer_meta_data);
    // writer2.write(&p_buffer[0]);

    // std::cout << my_buffer[0] << " " << my_buffer[1] << " " << my_buffer[2] << " " << my_buffer[3] << std::endl;

    // delete size;
    // return 0;

    struct jpeg_compress_struct cinfo;
    struct jpegutils::my_error_mgr jerr;

    // Initialization.

    cinfo.err = jpeg_std_error(&jerr.pub);
        jerr.pub.output_message = my_output_message;
        jerr.pub.error_exit = my_error_exit;
        if(setjmp(jerr.setjmp_buffer))
        {
           std::cout << "Error! Unable to intialize error handler." << std::endl;
        }

        jpeg_create_compress(&cinfo);

        cinfo.image_width  = (int) writer_meta_data.size[1];
        cinfo.image_height = (int) writer_meta_data.size[0]; 

        cinfo.in_color_space = writer_meta_data.imageType;
        if ( writer_meta_data.imageType == JCS_RGB ) {
            cinfo.input_components = 3;
        } else {
            // JCS_GRAYSCALE
            cinfo.input_components = 1;
        }


        /* Set default compression parameters. */
        jpeg_set_defaults(&cinfo);
        jpeg_set_quality(&cinfo, writer_meta_data.quality, TRUE);
        
        unsigned char * outbuffer;
        unsigned long outsize = 0;

        jpeg_mem_dest (&cinfo, &outbuffer, &outsize);

        if (writer_meta_data.compMode == MEX_JPEG_LOSSLESS) {
            jpeg_simple_lossless(&cinfo, 1, 0);
        }

        jpeg_start_compress(&cinfo, TRUE);

        JSAMPARRAY buffer;
    
    size_t height = size[0];
    size_t width  = size[1];
    
    JDIMENSION row_stride = boost::numeric_cast<JDIMENSION>(width * 3);    
    buffer = (cinfo.mem->alloc_sarray)
        ((j_common_ptr) &cinfo, JPOOL_IMAGE, row_stride, 1);
    
    /* jpeg_write_scanlines expects an array of pointers to scanlines.
     * Here the array is only one element long, but you could pass
     * more than one scanline at a time if that's more convenient.
     */
    
    while (cinfo.next_scanline < cinfo.image_height) {
        /* construct a buffer which contains the next scanline with data in
         * RGBRGBRGB... format */
        JDIMENSION i = cinfo.next_scanline;
        for(JDIMENSION j = 0; j < cinfo.image_width; j++)
        {
            buffer[0][3*j]   = uint8Data[i + (j*height)]; /* Red */
            buffer[0][3*j+1] = uint8Data[i + (j*height) + (height * width)]; /* Green */
            buffer[0][3*j+2] = uint8Data[i + (j*height) + (2 * height * width)]; /* Blue */
        }
        (void) jpeg_write_scanlines(&cinfo, buffer, 1);
    } 

    std::cout << outsize << std::endl;
    for(int i = 0; i < 4096; i++) {
        std::cout << static_cast<int>(outbuffer[i]) << " ";
        if(i % 50 == 0) {
            std::cout << std::endl;
        }
    }

    jpeg_finish_compress(&cinfo);
    jpeg_destroy_compress(&cinfo);

    // std::cout << std::endl << std::endl;

    // std::string infilename2 = "Resources/rosew1.jpg";
    
    // FILE * pFile2 = std::fopen(infilename.c_str(), "rb");
    // if (pFile2 == NULL) {
    //     fputs("File error", stderr);
    //     exit(1);
    // }

    // // obtain file size:
    // fseek(pFile2, 0, SEEK_END);
    // size_t lSize2 = ftell(pFile);
    // rewind(pFile2);

    // // allocate memory to contain the whole file:
    // if (p_buffer.empty()) {
    //     fputs("Memory error", stderr);
    //     exit(2);
    // }

    // std::vector<unsigned char> p_buffer2(lSize2);

    // // copy the file contents into the buffer:
    // size_t result2 = fread(&p_buffer2[0], sizeof(unsigned char), lSize2, pFile2);
    // if (result != lSize2) {
    //     fputs("Reading error", stderr);
    //     exit(3);
    // }

    // fclose(pFile2);

    // for(int i = 0; i < lSize2; i++) {
    //     std::cout << static_cast<int>(p_buffer2[i]) << " ";
    //     if(i % 50 == 0) {
    //         std::cout << std::endl;
    //     }
    // }

    return 0;
}