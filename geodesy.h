/*
 * Библиотека перевода геодезических координат
 * 
 * История версий:
 * 
 * Версия 3.0 (текущая), 03.02.2021
 *     Удалены функции ecef2ecef_xxx в пользу констант преобразования
 *
 *     Для функции ecef2ecef добавлены константы преобразования:
 *         wgs84_to_pz9011
 *         pz9011_to_wgs84
 *         pz9011_to_pz90
 *         pz90_to_pz9011
 *         sk42_to_pz9011
 *         pz9011_to_sk42
 *
 * Версия 2.0, 17.09.2020
 *     Введены doxygen-комментарии
 *
 *     Добавлены функции перевода:
 *         geodetic2gk_zone        - определение предпочтительной зоны в проекции Гаусса-Крюгера
 *         gk2geodetic             - перевод из проекции Гаусса-Крюгера в геодезическую СК (СК-42)
 *         ecef2ecef_wgs84_pz90    - перевод геоцентрических координат из WGS-84 в ПЗ-90
 *         ecef2ecef_ck42_pz90     - перевод геоцентрических координат из СК-42  в ПЗ-90
 *         ecef2ecef_ck42_wgs84    - перевод геоцентрических координат из СК-42  в WGS-84
 *         ecef2ecef_wgs84_ck42    - перевод геоцентрических координат из WGS-84 в СК-42
 *         ecef2ecef_pz90_ck42     - перевод геоцентрических координат из ПЗ-90  в СК-42
 *         ecef2ecef_pz90_wgs84    - перевод геоцентрических координат из ПЗ-90  в WGS-84
 *
 * Версия 1.0, 11.09.2020
 *     Добавлены основные функции перевода:
 *         geodetic2ecef           - перевод из геодезической СК в геоцентрическую с указанием параметров эллипсоида
 *         ecef2geodetic           - перевод из геоцентрической СК в геодезическую с указанием параметров эллипсоида
 *         ecef2enu                - перевод из геоцентрической СК в топоцентрическую с указанием точки стояния и эллипсоида
 *         enu2ecef                - перевод из топоцентрической СК в геоцентрическую с указанием точки стояния и эллипсоида
 *         geodetic2enu            - перевод из геодезической СК в топоцентрическую с указанием точки стояния и эллипсоида
 *         enu2geodetic            - перевод из топоцентрической СК в геодезическую с указанием точки стояния и эллипсоида
 *         geodetic2gk             - перевод из геодезической СК (СК-42) в проекцию Гаусса-Крюгера
 *         ecef2ecef               - перевод из геоцентрической СК в другую геоцентрическую СК
 */

#ifndef GEODESY_H
#define GEODESY_H

#include <stdint.h>
#include <stddef.h>
#include <math.h>

#ifdef __cplusplus
extern "C" {
#endif

typedef struct
{
    double a; /**< Большая полуось, метры */
    double b; /**< Малая полуось, метры */
    double f; /**< Коэффициент сжатия: 1/f = a / (a - b) */
} ellipsoid_t;

typedef struct
{
    double dx; /**< Линейный параметр трансформирования по х */
    double dy; /**< Линейный параметр трансформирования по y */
    double dz; /**< Линейный параметр трансформирования по z */
    double wx; /**< Угловой параметр трансформирования по х */
    double wy; /**< Угловой параметр трансформирования по y */
    double wz; /**< Угловой параметр трансформирования по z */
    double m;  /**< Масштабный параметр трансформирования */
} ecefconv_t;

extern const ellipsoid_t sk42;           /**< Параметры эллипсоида СК-42 */
extern const ellipsoid_t wgs84;          /**< Параметры эллипсоида WGS-84 */
extern const ellipsoid_t pz90;           /**< Параметры эллипсоида ПЗ-90 */
extern const ellipsoid_t pz9011;         /**< Параметры эллипсоида ПЗ-90.11 */

extern const ecefconv_t wgs84_to_pz9011; /**< Константы преобразования из wgs84 в ПЗ-90.11 */
extern const ecefconv_t pz9011_to_wgs84; /**< Константы преобразования из ПЗ-90.11 в wgs84 */
extern const ecefconv_t pz9011_to_pz90;  /**< Константы преобразования из ПЗ-90.11 в ПЗ-90 */
extern const ecefconv_t pz90_to_pz9011;  /**< Константы преобразования из ПЗ-90 в ПЗ-90.11 */
extern const ecefconv_t sk42_to_pz9011;  /**< Константы преобразования из СК-42 в ПЗ-90.11 */
extern const ecefconv_t pz9011_to_sk42;  /**< Константы преобразования из ПЗ-90.11 в СК-42 */

/**
 * @brief Перевод геодезических коодринат (ШДВ) в геоцентрические (XYZ)
 * @param [in]  lat Широта
 * @param [in]  lon Долгота
 * @param [in]  alt Высота над уровнем эллипсоида
 * @param [out] x   Рассчитанная координата X
 * @param [out] y   Рассчитанная координата Y
 * @param [out] z   Рассчитанная координата Z
 * @param [in]  ell Параметры эллипсоида
 */
void geodetic2ecef(double lat, double lon, double alt, double *x, double *y, double *z, const ellipsoid_t *ell);

/**
 * @brief Перевод геоцентрических коодринат (XYZ) в геодезические (ШДВ)
 * @param [in]  x   Координата X
 * @param [in]  y   Координата Y
 * @param [in]  z   Координата Z
 * @param [out] lat Рассчитанная широта
 * @param [out] lon Рассчитанная долгота
 * @param [out] alt Рассчитанная высота над эллипсоидом
 * @param [in]  ell Параметры эллипсоида
 */
void ecef2geodetic(double x, double y, double z, double *lat, double *lon, double *alt, const ellipsoid_t *ell);

/**
 * @brief Перевод геоцентрических коодринат (XYZ) в топоцентрические (ENU)
 * @param [in]  x     Координата X
 * @param [in]  y     Координата Y
 * @param [in]  z     Координата Z
 * @param [in]  lat0  Широта начала координат топоцентрической системы
 * @param [in]  lon0  Долгота начала координат топоцентрической системы
 * @param [in]  h0    Высота начала координат топоцентрической системы
 * @param [out] east  Рассчитанная координата в топоцентрической системе в направлении востока
 * @param [out] north Рассчитанная координата в топоцентрической системе в направлении севера
 * @param [out] up    Рассчитанная координата в топоцентрической системе в направлении вверх
 * @param [in]  ell   Параметры эллипсоида
 */
void ecef2enu(double x, double y, double z, double lat0, double lon0, double h0, double *east, double *north, double *up, const ellipsoid_t *ell);

/**
 * @brief Перевод топоцентрических координат (ENU) в геоцентрические (XYZ)
 * @param [in]  e1   Координата в топоцентрической системе в направлении востока
 * @param [in]  n1   Координата в топоцентрической системе в направлении севера
 * @param [in]  u1   Координата в топоцентрической системе в направлении вверх
 * @param [in]  lat0 Широта начала координат топоцентрической системы
 * @param [in]  lon0 Долгота начала координат топоцентрической системы
 * @param [in]  h0   Высота начала координат топоцентрической системы
 * @param [out] x    Рассчитанная координата X
 * @param [out] y    Рассчитанная координата Y
 * @param [out] z    Рассчитанная координата Z
 * @param [in]  ell  Параметры эллипсоида
 */
void enu2ecef(double e1, double n1, double u1, double lat0, double lon0, double h0, double *x, double *y, double *z, const ellipsoid_t *ell);

/**
 * @brief Перевод геодезических координат (ШДВ) в топоцентрические (ENU)
 * @param [in]  lat   Широта
 * @param [in]  lon   Долгота
 * @param [in]  h     Высота над уровнем эллипсоида
 * @param [in]  lat0  Широта начала координат топоцентрической системы
 * @param [in]  lon0  Долгота начала координат топоцентрической системы
 * @param [in]  h0    Высота начала координат топоцентрической системы
 * @param [out] east  Рассчитанная координата в топоцентрической системе в направлении востока
 * @param [out] north Рассчитанная координата в топоцентрической системе в направлении севера
 * @param [out] up    Рассчитанная координата в топоцентрической системе в направлении вверх
 * @param [in]  ell   Параметры эллипсоида
 */
void geodetic2enu(double lat, double lon, double h, double lat0, double lon0, double h0, double *east, double *north, double *up, const ellipsoid_t *ell);

/**
 * @brief Перевод топоцентрических (ENU) координат в геодезические (ШДВ)
 * @param [in]  east  Координата в топоцентрической системе в направлении востока
 * @param [in]  north Координата в топоцентрической системе в направлении севера
 * @param [in]  up    Координата в топоцентрической системе в направлении вверх
 * @param [in]  lat0  Широта начала координат топоцентрической системы
 * @param [in]  lon0  Долгота начала координат топоцентрической системы
 * @param [in]  h0    Высота начала координат топоцентрической системы
 * @param [out] lat   Рассчитанная широта
 * @param [out] lon   Рассчитанная долгота
 * @param [out] h     Рассчитанная высота над эллипсоидом
 * @param [in]  ell   Параметры эллипсоида
 */
void enu2geodetic(double east, double north, double up, double lat0, double lon0, double h0, double *lat, double *lon, double *h, const ellipsoid_t *ell);

/**
 * @brief Перевод из геодезической СК в проекцию Гаусса - Крюгера
 * @param [in]  lat  Широта
 * @param [in]  lon  Долгота
 * @param [in]  zone Номер зоны для вычислений проекции
 * @param [out] x_gk Рассчитанная координата x
 * @param [out] y_gk Рассчитанная координата y
 */
void geodetic2gk(double lat, double lon, int zone, double *x_gk, double *y_gk);

/**
 * @brief Возвращает номер предпочтительной для использования зоны в проекции Гаусса-Крюгера
 * @param [in] lon Долгота [0, 2pi)
 * @return Номер зоны в проекции Гаусса-Крюгера
 */
int geodetic2gk_zone(double lon);

/**
 * @brief Перевод из проекции Гаусса - Крюгера в геодезическую СК
 * @param [in]  x    Абсцисса
 * @param [in]  y    Ордината
 * @param [out] lat  Рассчитанная широта x
 * @param [out] lon  Рассчитанная долгота y
 */
void gk2geodetic(double x, double y, double* lat, double* lon);

/**
 * @brief Перевод из проекции Гаусса - Крюгера в геодезическую СК
 * @param [in]  x    Абсцисса
 * @param [in]  y    Ордината
 * @param [in]  zone Номер зоны для вычисленной проекции
 * @param [out] lat  Рассчитанная широта x
 * @param [out] lon  Рассчитанная долгота y
 */
void gkz2geodetic(double x, double y, int zone, double* lat, double* lon);

/**
 * @brief Возвращает номер зоны в проекции Гаусса-Крюгера
 * @param [in] y Ордината [0, 60 * 1e6)
 * @return Номер зоны в проекции Гаусса-Крюгера
 */
int gk_zone(double y);

/**
 * @brief Перевод из одной геоцентрической системы координат в другую
 * @param [in]  src_x координата x в изначальной системе координат
 * @param [in]  src_y координата y в изначальной системе координат
 * @param [in]  src_z координата z в изначальной системе координат
 * @param [out] dst_x координата x в новой системе координат
 * @param [out] dst_y координата y в новой системе координат
 * @param [out] dst_z координата z в новой системе координат
 * @param [in] conv Константы преобразования
 */
void ecef2ecef(double src_x, double src_y, double src_z, double *dst_x, double *dst_y, double *dst_z, const ecefconv_t *conv);

#ifdef __cplusplus
}
#endif

#endif /* GEODESY_H */
