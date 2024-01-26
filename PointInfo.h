#ifndef POINT_INFO_H_
#define POINT_INFO_H_

#include <cstdint>
/* Определеим Структуры и функции для управления информацией о точке в расчетной области */

namespace PointInfo2D
{
    struct Info
    {
    };

    struct BoundInfo
    {
    };

    struct AreaInfo
    {
    };

    struct Comand
    {
    };

    /* Функции управления битовыми структурами */
    /* Private Section */

    /* Очистка структуры */
    void ClearInfo(Info &info);

    /* Установка фиктивный/не фиктивный узел
        @param
        Info& info - информационная структура
        uint8_t val true - не фиктивный false - фиктивный
        @return void
    */
    void SetFictitiousBit(Info &info, uint8_t val);

    /* Установка Количества различныйх типов областей которым принадлежит точка
        @param
        Info& info - информационная структура
        uint8_t val - диапазон [0;8]
        @return void
    */
    void SetAreaCountBits(Info &info, uint32_t val);

    /* Установка Количества различных типов границ к которым примыкает данная точка
        @param
        Info& info - информационная структура
        uint8_t val - диапазон [0;4]
        @return void
    */
    void SetBoundCountBits(Info &info, uint32_t val);

    /* Получить количество различных областей к которым примыкает точка
        @param
        Info& info - информационная структура
        @return uint32_t - информациоя о том к скольки областям примыкает  диапазон [0;8]
        0 - ни к одной
        1 - к одной
        ...
    */
    uint32_t GetAreaCount(const Info &info);

    /* Получить количество различных областей к которым примыкает граница
        @param
        Info& info - информационная структура
        @return uint32_t - информациоя о том к скольки Границам примыкает(различным)  диапазон [0;4]
        0 - ни к какой
        1 - к одной
        2 -
        3 -
        4 -
    */
    uint32_t GetBoundCount(const Info &info);
    /*********************************************/

    /* Public section */

    /*
        @param
        Info& info
        uint32_t val - набор формул задающий значения на границе диапазон [0,15]
        uint32_t boundType - тип Ку диапазон  [0,3]
    */
    void SetBoundInfo(Info &info, uint32_t val, uint32_t boundType);
    /*
        @param
        Info& info
        uint32_t val - набор формул задающий значения в области определяющей параметры ДУ  диапазон [0,15]
    */
    void SetAreaInfo(Info &info, uint32_t val);

    /*
        @param
        const Info& info
        @return bool - false - фиктивный true - не фиктивный
    */
    bool IsFiFictitious(const Info &info);

    /*
        @param
        const Info& info
        @return bool - false - не граница true - граница
    */
    bool IsBound(const Info &info);

    /*
        @param
        const Info& info
        @return BoundInfo -структура содержащая информацию о типах Ку и формул которые их задают
    */
    BoundInfo GetBoundInfo(const Info &info);

    /*
        @param
        const Info& info
        @return AreaInfo -структура содержащая информацию о типах областей которым принадлежит точка и как следствие это определяет парметры ДУ
    */
    AreaInfo GetAreaInfo(const Info &info);
    /********************************************/

    /* Debug functions */
    void PrintInfo(const Info &info);

    void PrintBoundInfo(const BoundInfo &Bound);

    void PrintAreaInfo(const AreaInfo &Area);

};

namespace PointInfo3D
{
    struct Info
    {
    };

    struct BoundInfo
    {
    };

    struct AreaInfo
    {
    };

    struct Comand
    {
    };
};

#endif
