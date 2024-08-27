/*
 *  author: Suhas Vittal
 *  date:   26 August 2024
 * */

#ifndef CGEN_SC_COORD_h
#define CGEN_SC_COORD_h

namespace cgen {

typedef uint16_t coord_t;
const uint16_t PZERO = static_cast<uint16_t>(0);

inline coord_t make_coord(uint8_t x, uint8_t y) { 
    return static_cast<uint16_t>(x) | (static_cast<uint16_t>(y) << 8);
}

inline uint8_t gx16(const coord_t& p) { return p & 0x00ff; }
inline uint8_t gy16(const coord_t& p) { return p & 0xff00; }
inline uint8_t gx8(const coord_t& p) { return static_cast<uint8_t>( gx16(p) ); }
inline uint8_t gy8(const coord_t& p) { return static_cast<uint8_t>( gy16(p) ); }

template <uint16_t DY>
inline uint16_t dim(const coord_t& p) { return gx16(p)*DY + gy16(p); }

inline coord_t c_add(coord_t a, coord_t b) {
    return make_coord( gx8(a)+gx8(b), gy8(a)+gy8(b) );
}

inline coord_t c_sub(coord_t a, coord_t b) {
    return make_coord( gx8(a)-gx8(b), gy8(a)-gy8(b) );
}

}   // cgen

#endif  // CGEN_SC_COORD_h
