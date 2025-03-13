#define _CRT_SECURE_NO_WARNINGS 1
#include <vector>

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"

#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"

#include <iostream>
#include <random>
#include <stdexcept>
#include <string>
#include <stdio.h>
#include <algorithm>
#include <list>

std::default_random_engine engine;
std::uniform_real_distribution<double> uniform(0, 1);

static inline double sqr(double x) { return x * x; }

class Vector
{
public:
    explicit Vector(double x = 0, double y = 0, double z = 0)
    {
        coord[0] = x;
        coord[1] = y;
        coord[2] = z;
    }
    double &operator[](int i) { return coord[i]; }
    double operator[](int i) const { return coord[i]; }
    Vector operator-() const { return Vector(-coord[0], -coord[1], -coord[2]); }

    Vector &operator+=(const Vector &v)
    {
        coord[0] += v[0];
        coord[1] += v[1];
        coord[2] += v[2];
        return *this;
    }

    double norm2() const
    {
        return sqr(coord[0]) + sqr(coord[1]) + sqr(coord[2]);
    }

    Vector normalize()
    {
        Vector v(coord[0], coord[1], coord[2]);
        coord[0] /= sqrt(v.norm2());
        coord[1] /= sqrt(v.norm2());
        coord[2] /= sqrt(v.norm2());
        return *this;
    }

    double coord[3];
};

Vector operator+(const Vector &a, const Vector &b)
{
    return Vector(a[0] + b[0], a[1] + b[1], a[2] + b[2]);
}
Vector operator-(const Vector &a, const Vector &b)
{
    return Vector(a[0] - b[0], a[1] - b[1], a[2] - b[2]);
}
Vector operator*(const Vector &a, double b)
{
    return Vector(a[0] * b, a[1] * b, a[2] * b);
}
Vector operator*(double a, const Vector &b)
{
    return Vector(a * b[0], a * b[1], a * b[2]);
}
Vector operator*(const Vector &a, const Vector &b)
{
    return Vector(a[0] * b[0], a[1] * b[1], a[2] * b[2]);
}
Vector operator/(const Vector &a, double b)
{
    return Vector(a[0] / b, a[1] / b, a[2] / b);
}
Vector operator+(const Vector &a, double b)
{
    return Vector(a[0] + b, a[1] + b, a[2] + b);
}

double dot(const Vector &a, const Vector &b)
{
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

Vector cross(const Vector &a, const Vector &b)
{
    return Vector(a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]);
}

Vector random_cos(const Vector &v)
{
    double random1 = uniform(engine);
    double random2 = uniform(engine);

    Vector direction_rebond(cos(2 * M_PI * random1) * sqrt(1 - random2), sin(2 * M_PI * random1) * sqrt(1 - random2), sqrt(random2));
    Vector randomV(uniform(engine) - 0.5, uniform(engine) - 0.5, uniform(engine) - 0.5);
    Vector tangeante1 = cross(v, randomV);
    tangeante1.normalize();
    Vector tangeante2 = cross(tangeante1, v);

    Vector direction_aleatoire = direction_rebond[2] * v + direction_rebond[0] * tangeante1 + direction_rebond[1] * tangeante2;
    return direction_aleatoire;
}

class Rayon
{
public:
    Vector origin;
    Vector u;

    Rayon()
    {
        origin = Vector(0.0, 0.0, 0.0);
        u = Vector(0.0, 0.0, 0.0);
    }

    Rayon(const Vector &o, const Vector &d)
    {
        origin = o;
        u = d;
    }
};

class Geometry
{
public:
    Vector albedo;
    bool mirroir;
    bool transparent;
    bool procedural;
    int index_procedural;

    Geometry()
    {
        albedo = Vector(1.0, 1.0, 1.0);
        mirroir = false;
        transparent = false;
        procedural = false;
        index_procedural = 0;
    }

    Geometry(const Vector &alb, bool m = false, bool t = false, bool proc = false, int index_proc = 0)
    {
        albedo = alb;
        mirroir = m;
        transparent = t;
        procedural = proc;
        index_procedural = index_proc;
    }

    virtual bool intersect(const Rayon &r, Vector &intersect, Vector &normale, double &t, Vector &albedo) const = 0;
};

class Sphere : public Geometry
{
public:
    float radius;
    Vector center;

    Sphere()
    {
        center = Vector(0.0, 0.0, 0.0);
        radius = (1.0);
        albedo = Vector(1.0, 1.0, 1.0);
        mirroir = false;
        transparent = false;
    };

    Sphere(const Vector &c, float r, const Vector &alb, bool m = false, bool t = false, bool proc = false, int index_proc = 0) : Geometry::Geometry(alb, m, t, proc, index_proc)
    {
        center = c;
        radius = r;
    };

    bool intersect(const Rayon &r, Vector &intersect, Vector &normale, double &t, Vector &albedo_sphere) const
    {
        Vector v = r.origin - center;
        double a = dot(r.u, r.u);
        double b = 2 * dot(v, r.u);
        double c = dot(v, v) - radius * radius;

        double delta = b * b - 4 * a * c;

        if (delta < 0)
        {
            return false;
        }

        double t1 = (-b + sqrt(delta)) / (2 * a);
        double t2 = (-b - sqrt(delta)) / (2 * a);
        if (t1 < 0 && t2 < 0)
        {
            return false;
        }
        if (t1 < 0)
        {
            t = t2;
        }
        else if (t2 < 0)
        {
            t = t1;
        }
        else
        {
            t = std::min(t1, t2);
        }

        intersect = r.origin + t * r.u;
        normale = (intersect - center) * (1. / radius);

        if (procedural)
        {
            // rendu parquet
            if (index_procedural == 0)
            {
                Vector value = Vector((intersect[0]), (intersect[1]), (intersect[2]));
                if ((int)(intersect[0]) >= 0)
                {
                    value[0] = value[0] + 1;
                }
                if ((int)(intersect[2]) >= 0)
                {
                    value[2] += 1;
                }
                if ((((int)(value[0]) + (int)(value[2])) % 3 == 0))
                {
                    albedo_sphere = Vector(0.43, .14, 0.05);
                }
                else
                {
                    albedo_sphere = Vector(.6, .3, .1);
                }
            }
            else
            {
                double noise = sin(10 * intersect[0]) * sin(10 * intersect[1]) * sin(10 * intersect[2]);
                albedo_sphere = Vector(0.5, 0.5, 0.9) * (1 + noise);
            }
        }
        else
            albedo_sphere = albedo;

        return true;
    }
};

class TriangleIndices
{
public:
    TriangleIndices(int vtxi = -1, int vtxj = -1, int vtxk = -1, int ni = -1, int nj = -1, int nk = -1, int uvi = -1, int uvj = -1, int uvk = -1, int group = -1, bool added = false) : vtxi(vtxi), vtxj(vtxj), vtxk(vtxk), uvi(uvi), uvj(uvj), uvk(uvk), ni(ni), nj(nj), nk(nk), group(group){};
    int vtxi, vtxj, vtxk; // indices within the vertex coordinates array
    int uvi, uvj, uvk;    // indices within the uv coordinates array
    int ni, nj, nk;       // indices within the normals array
    int group;            // face group
};

class BBox
{
public:
    Vector pMin, pMax;

    BBox(){};
    bool intersect(const Rayon &r) const
    {

        Vector invDir(1.0 / r.u[0], 1.0 / r.u[1], 1.0 / r.u[2]);

        double t1x = (pMin[0] - r.origin[0]) * invDir[0];
        double t2x = (pMax[0] - r.origin[0]) * invDir[0];

        double tinx = std::min(t1x, t2x);
        double toutx = std::max(t1x, t2x);

        double t1y = (pMin[1] - r.origin[1]) * invDir[1];
        double t2y = (pMax[1] - r.origin[1]) * invDir[1];

        double tiny = std::min(t1y, t2y);
        double touty = std::max(t1y, t2y);

        double t1z = (pMin[2] - r.origin[2]) * invDir[2];
        double t2z = (pMax[2] - r.origin[2]) * invDir[2];

        double tinz = std::min(t1z, t2z);
        double toutz = std::max(t1z, t2z);

        double tin = std::max(tinx, std::max(tiny, tinz));
        double tout = std::min(toutx, std::min(touty, toutz));

        if (tout < 0)
        {
            return false;
        }

        if (tout > tin)
        {
            return true;
        }

        return false;
    };
};

class BVH
{
public:
    int debut, fin;
    BBox bbox;

    BVH *fg;
    BVH *fd;
    BVH(const BBox &bbox, int index_init = -1, int index_end = -1, BVH *fg = NULL, BVH *fd = NULL) : bbox(bbox), fg(fg), fd(fd), debut(index_init), fin(index_end) {}
};

class TriangleMesh : public Geometry
{
public:
    ~TriangleMesh(){};
    TriangleMesh() : Geometry(){};

    BVH *root;

    void readOBJ(const char *obj)
    {

        char matfile[255];
        char grp[255];

        FILE *f;
        f = fopen(obj, "r");
        int curGroup = -1;
        while (!feof(f))
        {
            char line[255];
            if (!fgets(line, 255, f))
                break;

            std::string linetrim(line);
            linetrim.erase(linetrim.find_last_not_of(" \r\t") + 1);
            strcpy(line, linetrim.c_str());

            if (line[0] == 'u' && line[1] == 's')
            {
                sscanf(line, "usemtl %[^\n]\n", grp);
                curGroup++;
            }

            if (line[0] == 'v' && line[1] == ' ')
            {
                Vector vec;

                Vector col;
                if (sscanf(line, "v %lf %lf %lf %lf %lf %lf\n", &vec[0], &vec[1], &vec[2], &col[0], &col[1], &col[2]) == 6)
                {
                    col[0] = std::min(1., std::max(0., col[0]));
                    col[1] = std::min(1., std::max(0., col[1]));
                    col[2] = std::min(1., std::max(0., col[2]));

                    vertices.push_back(vec);
                    vertexcolors.push_back(col);
                }
                else
                {
                    sscanf(line, "v %lf %lf %lf\n", &vec[0], &vec[1], &vec[2]);
                    vertices.push_back(vec);
                }
            }
            if (line[0] == 'v' && line[1] == 'n')
            {
                Vector vec;
                sscanf(line, "vn %lf %lf %lf\n", &vec[0], &vec[1], &vec[2]);
                normals.push_back(vec);
            }
            if (line[0] == 'v' && line[1] == 't')
            {
                Vector vec;
                sscanf(line, "vt %lf %lf\n", &vec[0], &vec[1]);
                uvs.push_back(vec);
            }
            if (line[0] == 'f')
            {
                TriangleIndices t;
                int i0, i1, i2, i3;
                int j0, j1, j2, j3;
                int k0, k1, k2, k3;
                int nn;
                t.group = curGroup;

                char *consumedline = line + 1;
                int offset;

                nn = sscanf(consumedline, "%u/%u/%u %u/%u/%u %u/%u/%u%n", &i0, &j0, &k0, &i1, &j1, &k1, &i2, &j2, &k2, &offset);
                if (nn == 9)
                {
                    if (i0 < 0)
                        t.vtxi = vertices.size() + i0;
                    else
                        t.vtxi = i0 - 1;
                    if (i1 < 0)
                        t.vtxj = vertices.size() + i1;
                    else
                        t.vtxj = i1 - 1;
                    if (i2 < 0)
                        t.vtxk = vertices.size() + i2;
                    else
                        t.vtxk = i2 - 1;
                    if (j0 < 0)
                        t.uvi = uvs.size() + j0;
                    else
                        t.uvi = j0 - 1;
                    if (j1 < 0)
                        t.uvj = uvs.size() + j1;
                    else
                        t.uvj = j1 - 1;
                    if (j2 < 0)
                        t.uvk = uvs.size() + j2;
                    else
                        t.uvk = j2 - 1;
                    if (k0 < 0)
                        t.ni = normals.size() + k0;
                    else
                        t.ni = k0 - 1;
                    if (k1 < 0)
                        t.nj = normals.size() + k1;
                    else
                        t.nj = k1 - 1;
                    if (k2 < 0)
                        t.nk = normals.size() + k2;
                    else
                        t.nk = k2 - 1;
                    indices.push_back(t);
                }
                else
                {
                    nn = sscanf(consumedline, "%u/%u %u/%u %u/%u%n", &i0, &j0, &i1, &j1, &i2, &j2, &offset);
                    if (nn == 6)
                    {
                        if (i0 < 0)
                            t.vtxi = vertices.size() + i0;
                        else
                            t.vtxi = i0 - 1;
                        if (i1 < 0)
                            t.vtxj = vertices.size() + i1;
                        else
                            t.vtxj = i1 - 1;
                        if (i2 < 0)
                            t.vtxk = vertices.size() + i2;
                        else
                            t.vtxk = i2 - 1;
                        if (j0 < 0)
                            t.uvi = uvs.size() + j0;
                        else
                            t.uvi = j0 - 1;
                        if (j1 < 0)
                            t.uvj = uvs.size() + j1;
                        else
                            t.uvj = j1 - 1;
                        if (j2 < 0)
                            t.uvk = uvs.size() + j2;
                        else
                            t.uvk = j2 - 1;
                        indices.push_back(t);
                    }
                    else
                    {
                        nn = sscanf(consumedline, "%u %u %u%n", &i0, &i1, &i2, &offset);
                        if (nn == 3)
                        {
                            if (i0 < 0)
                                t.vtxi = vertices.size() + i0;
                            else
                                t.vtxi = i0 - 1;
                            if (i1 < 0)
                                t.vtxj = vertices.size() + i1;
                            else
                                t.vtxj = i1 - 1;
                            if (i2 < 0)
                                t.vtxk = vertices.size() + i2;
                            else
                                t.vtxk = i2 - 1;
                            indices.push_back(t);
                        }
                        else
                        {
                            nn = sscanf(consumedline, "%u//%u %u//%u %u//%u%n", &i0, &k0, &i1, &k1, &i2, &k2, &offset);
                            if (i0 < 0)
                                t.vtxi = vertices.size() + i0;
                            else
                                t.vtxi = i0 - 1;
                            if (i1 < 0)
                                t.vtxj = vertices.size() + i1;
                            else
                                t.vtxj = i1 - 1;
                            if (i2 < 0)
                                t.vtxk = vertices.size() + i2;
                            else
                                t.vtxk = i2 - 1;
                            if (k0 < 0)
                                t.ni = normals.size() + k0;
                            else
                                t.ni = k0 - 1;
                            if (k1 < 0)
                                t.nj = normals.size() + k1;
                            else
                                t.nj = k1 - 1;
                            if (k2 < 0)
                                t.nk = normals.size() + k2;
                            else
                                t.nk = k2 - 1;
                            indices.push_back(t);
                        }
                    }
                }

                consumedline = consumedline + offset;

                while (true)
                {
                    if (consumedline[0] == '\n')
                        break;
                    if (consumedline[0] == '\0')
                        break;
                    nn = sscanf(consumedline, "%u/%u/%u%n", &i3, &j3, &k3, &offset);
                    TriangleIndices t2;
                    t2.group = curGroup;
                    if (nn == 3)
                    {
                        if (i0 < 0)
                            t2.vtxi = vertices.size() + i0;
                        else
                            t2.vtxi = i0 - 1;
                        if (i2 < 0)
                            t2.vtxj = vertices.size() + i2;
                        else
                            t2.vtxj = i2 - 1;
                        if (i3 < 0)
                            t2.vtxk = vertices.size() + i3;
                        else
                            t2.vtxk = i3 - 1;
                        if (j0 < 0)
                            t2.uvi = uvs.size() + j0;
                        else
                            t2.uvi = j0 - 1;
                        if (j2 < 0)
                            t2.uvj = uvs.size() + j2;
                        else
                            t2.uvj = j2 - 1;
                        if (j3 < 0)
                            t2.uvk = uvs.size() + j3;
                        else
                            t2.uvk = j3 - 1;
                        if (k0 < 0)
                            t2.ni = normals.size() + k0;
                        else
                            t2.ni = k0 - 1;
                        if (k2 < 0)
                            t2.nj = normals.size() + k2;
                        else
                            t2.nj = k2 - 1;
                        if (k3 < 0)
                            t2.nk = normals.size() + k3;
                        else
                            t2.nk = k3 - 1;
                        indices.push_back(t2);
                        consumedline = consumedline + offset;
                        i2 = i3;
                        j2 = j3;
                        k2 = k3;
                    }
                    else
                    {
                        nn = sscanf(consumedline, "%u/%u%n", &i3, &j3, &offset);
                        if (nn == 2)
                        {
                            if (i0 < 0)
                                t2.vtxi = vertices.size() + i0;
                            else
                                t2.vtxi = i0 - 1;
                            if (i2 < 0)
                                t2.vtxj = vertices.size() + i2;
                            else
                                t2.vtxj = i2 - 1;
                            if (i3 < 0)
                                t2.vtxk = vertices.size() + i3;
                            else
                                t2.vtxk = i3 - 1;
                            if (j0 < 0)
                                t2.uvi = uvs.size() + j0;
                            else
                                t2.uvi = j0 - 1;
                            if (j2 < 0)
                                t2.uvj = uvs.size() + j2;
                            else
                                t2.uvj = j2 - 1;
                            if (j3 < 0)
                                t2.uvk = uvs.size() + j3;
                            else
                                t2.uvk = j3 - 1;
                            consumedline = consumedline + offset;
                            i2 = i3;
                            j2 = j3;
                            indices.push_back(t2);
                        }
                        else
                        {
                            nn = sscanf(consumedline, "%u//%u%n", &i3, &k3, &offset);
                            if (nn == 2)
                            {
                                if (i0 < 0)
                                    t2.vtxi = vertices.size() + i0;
                                else
                                    t2.vtxi = i0 - 1;
                                if (i2 < 0)
                                    t2.vtxj = vertices.size() + i2;
                                else
                                    t2.vtxj = i2 - 1;
                                if (i3 < 0)
                                    t2.vtxk = vertices.size() + i3;
                                else
                                    t2.vtxk = i3 - 1;
                                if (k0 < 0)
                                    t2.ni = normals.size() + k0;
                                else
                                    t2.ni = k0 - 1;
                                if (k2 < 0)
                                    t2.nj = normals.size() + k2;
                                else
                                    t2.nj = k2 - 1;
                                if (k3 < 0)
                                    t2.nk = normals.size() + k3;
                                else
                                    t2.nk = k3 - 1;
                                consumedline = consumedline + offset;
                                i2 = i3;
                                k2 = k3;
                                indices.push_back(t2);
                            }
                            else
                            {
                                nn = sscanf(consumedline, "%u%n", &i3, &offset);
                                if (nn == 1)
                                {
                                    if (i0 < 0)
                                        t2.vtxi = vertices.size() + i0;
                                    else
                                        t2.vtxi = i0 - 1;
                                    if (i2 < 0)
                                        t2.vtxj = vertices.size() + i2;
                                    else
                                        t2.vtxj = i2 - 1;
                                    if (i3 < 0)
                                        t2.vtxk = vertices.size() + i3;
                                    else
                                        t2.vtxk = i3 - 1;
                                    consumedline = consumedline + offset;
                                    i2 = i3;
                                    indices.push_back(t2);
                                }
                                else
                                {
                                    consumedline = consumedline + 1;
                                }
                            }
                        }
                    }
                }
            }
        }
        fclose(f);
    }

    BBox build_bbox(int triangle_init, int triangle_end)
    {
        BBox bbox;
        bbox.pMin = Vector(1E10, 1E10, 1E10);
        bbox.pMax = Vector(-1E10, -1E10, -1E10);

        for (int i = triangle_init; i < triangle_end; i++)
        {
            for (int j = 0; j < 3; j++)
            {
                bbox.pMin[j] = std::min(bbox.pMin[j], vertices[indices[i].vtxi][j]);
                bbox.pMin[j] = std::min(bbox.pMin[j], vertices[indices[i].vtxj][j]);
                bbox.pMin[j] = std::min(bbox.pMin[j], vertices[indices[i].vtxk][j]);

                bbox.pMax[j] = std::max(bbox.pMax[j], vertices[indices[i].vtxi][j]);
                bbox.pMax[j] = std::max(bbox.pMax[j], vertices[indices[i].vtxj][j]);
                bbox.pMax[j] = std::max(bbox.pMax[j], vertices[indices[i].vtxk][j]);
            }
        }
        return bbox;
    }

    void build_bvh(BVH *node, int debut, int fin)
    {
        node->bbox = build_bbox(debut, fin);
        node->debut = debut;
        node->fin = fin;
        node->fg = NULL;
        node->fd = NULL;

        Vector size = node->bbox.pMax - node->bbox.pMin;
        Vector mid = (node->bbox.pMin + node->bbox.pMax) / 2;

        int dimension;
        if (size[0] > size[1] && size[0] > size[2])
            dimension = 0;
        else if (size[1] > size[2])
            dimension = 1;
        else
            dimension = 2;

        int pivot = debut;
        for (int i = debut; i < fin; i++)
        {
            double mid_triangle = (vertices[indices[i].vtxi][dimension] + vertices[indices[i].vtxj][dimension] + vertices[indices[i].vtxk][dimension]) / 3;
            if (mid_triangle < mid[dimension])
            {
                std::swap(indices[i], indices[pivot]);
                pivot++;
            }
        }

        // On tolere 5 triangles par noeud
        if (fin - debut < 5 || pivot - debut < 2 || fin - pivot < 2)
            return;

        node->fg = new BVH(node->bbox);
        build_bvh(node->fg, debut, pivot);

        node->fd = new BVH(node->bbox);
        build_bvh(node->fd, pivot, fin);
    }

    bool intersect(const Rayon &r, Vector &intersect, Vector &normal, double &t, Vector &albedo_mesh) const
    {
        double new_t = 1E10;
        bool intersected = false;
        double best_alpha, best_beta;
        int best_index;

        if (!root->bbox.intersect(r))
            return false;

        std::list<BVH *> nodes;
        nodes.push_back(root);

        while (!nodes.empty())
        {
            const BVH *actual = nodes.back();
            nodes.pop_back();

            if (actual->fg)
            {
                if (actual->fg->bbox.intersect(r))
                {
                    nodes.push_back(actual->fg);
                }
                if (actual->fd->bbox.intersect(r))
                {
                    nodes.push_back(actual->fd);
                }
            }
            else
            {
                for (int i = actual->debut; i < actual->fin; i++)
                {
                    const Vector &A = vertices[indices[i].vtxi];
                    const Vector &B = vertices[indices[i].vtxj];
                    const Vector &C = vertices[indices[i].vtxk];

                    Vector e1 = B - A;
                    Vector e2 = C - A;
                    Vector localN = cross(e1, e2);
                    Vector OA_u = cross(A - r.origin, r.u);

                    double uN = dot(r.u, localN);
                    double beta = dot(e2, OA_u) / uN;
                    double gamma = -dot(e1, OA_u) / uN;
                    double t_local = dot(A - r.origin, localN) / uN;
                    double alpha = 1 - beta - gamma;

                    if (beta >= 0 && gamma >= 0 && beta <= 1 && gamma <= 1 && alpha >= 0 && t_local > 0)
                    {
                        intersected = true;
                        if (t_local < new_t)
                        {
                            new_t = t_local;
                            t = t_local;
                            best_index = i;
                            best_alpha = alpha;
                            best_beta = beta;

                            intersect = r.origin + t * r.u;
                            normal = alpha * normals[indices[i].ni] + beta * normals[indices[i].nj] + gamma * normals[indices[i].nk];
                        }
                    }
                }
            }
        }
        if (intersected)
        {
            if (indices[best_index].group > texture.size() || texture.size() == 0)
            {
                albedo_mesh = albedo;
            }
            else
            {
                Vector Uvp = uvs[indices[best_index].uvi] * best_alpha + uvs[indices[best_index].uvj] * best_beta + uvs[indices[best_index].uvk] * (1 - best_alpha - best_beta);
                Uvp[0] = fabs(Uvp[0]);
                Uvp[0] = Uvp[0] - floor(Uvp[0]);
                Uvp[1] = fabs(Uvp[1]);
                Uvp[1] = Uvp[1] - floor(Uvp[1]);
                Uvp[1] = 1 - Uvp[1];

                Uvp = Uvp * Vector(texture_width[indices[best_index].group], texture_height[indices[best_index].group], 0);

                int u = std::min((int)Uvp[0], texture_width[indices[best_index].group] - 1);
                int v = std::min((int)Uvp[1], texture_height[indices[best_index].group] - 1);
                int pixel = v * texture_width[indices[best_index].group] + u;

                albedo_mesh = Vector(texture[indices[best_index].group][3 * pixel], texture[indices[best_index].group][3 * pixel + 1], texture[indices[best_index].group][3 * pixel + 2]);
            }
        }
        return intersected;
    }

    void transform(double scale_factor, Vector translation, double rotation = 0, const Vector &axis = Vector(0, 0, 0))
    {
        for (int i = 0; i < vertices.size(); i++)
        {
            vertices[i] = vertices[i] * scale_factor + translation;
        }

        // Matrice de rotation
        double s = sin(rotation);
        double c = cos(rotation);
        for (int i = 0; i < vertices.size(); i++)
        {
            double x = vertices[i][0];
            double y = vertices[i][1];
            double z = vertices[i][2];
            vertices[i][0] = x * (c + sqr(axis[0]) * (1 - c)) + y * (axis[0] * axis[1] * (1 - c) - axis[2] * s) + z * (axis[0] * axis[2] * (1 - c) + axis[1] * s);
            vertices[i][1] = x * (axis[1] * axis[0] * (1 - c) + axis[2] * s) + y * (c + sqr(axis[1]) * (1 - c)) + z * (axis[1] * axis[2] * (1 - c) - axis[0] * s);
            vertices[i][2] = x * (axis[2] * axis[0] * (1 - c) - axis[1] * s) + y * (axis[2] * axis[1] * (1 - c) + axis[0] * s) + z * (c + sqr(axis[2]) * (1 - c));
        }
    }

    void load_texture(const char *file)
    {
        int w, h, n;
        texture.push_back(stbi_loadf(file, &w, &h, &n, 3));
        texture_width.push_back(w);
        texture_height.push_back(h);
    }

    std::vector<float *> texture;
    std::vector<int> texture_width, texture_height;
    std::vector<TriangleIndices> indices;
    std::vector<Vector> vertices;
    std::vector<Vector> normals;
    std::vector<Vector> uvs;
    std::vector<Vector> vertexcolors;
};

class Scene
{
public:
    std::vector<const Geometry *> objects;

    Scene() {}

    void addGeometry(const Geometry *s)
    {
        objects.push_back(s);
    }

    bool intersect(const Rayon &r, Vector &intersect, Vector &normale, int &objectId, double &t, Vector &albedo) const
    {
        t = 1E11;
        bool intersected = false;
        for (int i = 0; i < objects.size(); i++)
        {
            Vector inter, norm, local_albedo;
            double t_obj;
            if (objects[i]->intersect(r, inter, norm, t_obj, local_albedo) && t_obj < t)
            {
                t = t_obj;
                intersect = inter;
                normale = norm;
                objectId = i;
                intersected = true;
                albedo = local_albedo;
            }
        }
        return intersected;
    }

    Vector getColor(const Rayon &r, int rebond, bool wasDiffuse = false)
    {
        double I = 1E10;
        double EPSILON = 0.001;
        Vector intersection, normale, albedo;
        int objectId;
        double t;
        Vector color;
        if (rebond == 0)
        {
            color = Vector(0., 0., 0.);
        }
        else if (intersect(r, intersection, normale, objectId, t, albedo))
        {
            if (objectId == 0)
            {
                if (wasDiffuse)
                {
                    return Vector(0., 0., 0.);
                }
                else
                {
                    // Pour avoir une lumière de couleur
                    return objects[0]->albedo * I / (4. * sqr(M_PI * dynamic_cast<const Sphere *>(objects[0])->radius));
                }
            }
            else if (objects[objectId]->mirroir)
            {
                Vector reflected = r.u - 2 * dot(r.u, normale) * normale;
                Rayon r2(intersection + EPSILON * normale, reflected);

                color = getColor(r2, rebond - 1);
            }
            else if (objects[objectId]->transparent)
            {
                double n1 = 1.0;
                double n2 = 1.333;
                Vector normale_transparence(normale);

                if (dot(r.u, normale) > 0)
                {
                    n1 = 1.333;
                    n2 = 1.0;
                    normale_transparence = -normale;
                }

                double calcul_racine = 1 - sqr(n1 / n2) * (1 - sqr(dot(normale_transparence, r.u)));
                if (calcul_racine > 0)
                {
                    // Sans Fresnel
                    Vector reflected = (n1 / n2) * (r.u - dot(r.u, normale_transparence) * normale_transparence) - normale_transparence * sqrt(calcul_racine);
                    Rayon r2(intersection - EPSILON * normale_transparence, reflected);
                    color = getColor(r2, rebond - 1);

                    // Avec Fresnel
                    double k0 = sqr((n1 - n2) / (n1 + n2));
                    double R = k0 + (1 - k0) * std::pow(1 - std::abs(dot(r.u, normale_transparence)), 5);
                    double T0 = 1 - R;
                    Rayon reflect = Rayon(intersection + EPSILON * normale_transparence, r.u - 2 * dot(r.u, normale_transparence) * normale_transparence);
                    color = T0 * color + R * getColor(reflect, rebond - 1);
                }
            }
            else
            {
                // Lumière directe
                Vector center = dynamic_cast<const Sphere *>(objects[0])->center;
                double radius = dynamic_cast<const Sphere *>(objects[0])->radius;

                Vector VecLumiere = center - intersection;
                VecLumiere.normalize();
                Vector Nprime = random_cos(-VecLumiere);
                Vector Pprime = Nprime * radius + center;
                Vector wi = Pprime - intersection;
                double d2 = wi.norm2();
                wi.normalize();

                double L_wi = I / (4. * sqr(M_PI * radius));
                Vector directe = albedo * L_wi * std::max(0., dot(normale, wi)) * std::max(0., dot(Nprime, -wi)) * (radius * radius) / (d2 * std::max(0., dot(Nprime, -VecLumiere)));

                // Lumière indirecte

                Vector direction_aleatoire = random_cos(normale);
                Rayon alea(intersection + EPSILON * normale, direction_aleatoire);

                Vector indirecte = albedo * getColor(alea, rebond - 1, true);

                // On regarde comment le pixel est eclairé

                double tShadow;
                Vector intersection2, normale2;

                Rayon r2(intersection + EPSILON * normale, wi);
                bool intersected = intersect(r2, intersection2, normale2, objectId, tShadow, albedo);
                double test = sqr(tShadow + 0.1);

                if (intersected && sqr(tShadow + 0.1) < d2)
                {
                    color = indirecte;
                }
                else
                {
                    color = directe + indirecte;
                }
                return color;
            }
        }
        return color;
    }
};

int main()
{
    // Paramètres rendu
    int W = 512;
    int H = 512;
    int numberRays = 64;
    int nRebonds = 5;

    // à régler pour avoir l'effet caméra ou non
    bool effetCamera = false;
    double distanceMiseAuPoint = 55;
    double ouverture = 1;

    Vector camera(0., 0., 55.);
    double fov = 60 * M_PI / 180;
    std::vector<unsigned char> image(W * H * 3, 0);

    // Scene
    Scene scene;

    // Load Mesh
    TriangleMesh mesh;
    mesh.readOBJ("./Models/cat.obj");
    mesh.load_texture("./Models/cat_texture.png");
    mesh.transform(0.4, Vector(5, -10, 10), 7 * M_PI / 4, Vector(0, 1, 0));
    mesh.root = new BVH(mesh.build_bbox(0, mesh.indices.size()));
    mesh.build_bvh(mesh.root, 0, mesh.indices.size());

    // Load Dinosaur
    TriangleMesh mesh2;
    mesh2.readOBJ("./Models/dinosaur.obj");
    mesh2.load_texture("./Models/dinosaur_texture.png");
    mesh2.transform(0.3, Vector(25, -10, 0), M_PI / 8, Vector(0, 1, 0));
    mesh2.root = new BVH(mesh2.build_bbox(0, mesh2.indices.size()));
    mesh2.build_bvh(mesh2.root, 0, mesh2.indices.size());

    // La première sphere est la source lumineuse
    scene.addGeometry(new Sphere(Vector(-10.0, 20.0, 40.0), 10., Vector(1., 1., 1.), false, false));

    // Ajout de sphères
    // scene.addGeometry(new Sphere(Vector(0.0, 15.0, -10.0), 5.0, Vector(1, 0., 1), true));
    // scene.addGeometry(new Sphere(Vector(-15.0, 5.0, 10.0), 5.0, Vector(0., 1, 1), false, true));
    // scene.addGeometry(new Sphere(Vector(-20.0, 25.0, 0.0), 5.0, Vector(1., 1., 1.), false, false));
    scene.addGeometry(new Sphere(Vector(0.0, 0.0, -1000.), 940., Vector(0., 1., 0.), false));
    scene.addGeometry(new Sphere(Vector(0.0, 1000.0, 0.0), 940., Vector(1., 0., 0.), false));
    scene.addGeometry(new Sphere(Vector(0.0, 0.0, 1000.), 940, Vector(1., 0., 0.5), false));
    scene.addGeometry(new Sphere(Vector(0.0, -1000.0, 0.0), 990, Vector(0., 0., 1.0), false, false, true));
    scene.addGeometry(new Sphere(Vector(-1000.0, 0.0, 0.0), 940, Vector(0., 1., 1.), false));
    scene.addGeometry(new Sphere(Vector(1000.0, 0.0, 0.0), 940, Vector(1., 1., 0.), false));

    // scene.addGeometry(&mesh);
    // scene.addGeometry(&mesh2);

#pragma omp parallel for schedule(dynamic, 1)
    for (int i = 0; i < H; i++)
    {
        for (int j = 0; j < W; j++)
        {
            Vector lumiere = Vector(0., 0., 0.);
            for (int k = 0; k < numberRays; k++)
            {
                // Box Muller - Anti aliasing
                double u1 = uniform(engine);
                double u2 = uniform(engine);

                double R = sqrt(-2 * log(u1));

                double dx = R * cos(2 * M_PI * u2) * .25;
                double dy = R * sin(2 * M_PI * u2) * .25;

                Vector u(j - W / 2. + 0.5 + dx, -i + W / 2. - 0.5 + dy, -W / (2. * tan(fov / 2.)));
                u.normalize();

                Rayon rayonCamera(camera, u);

                if (effetCamera)
                {
                    // Effet camera
                    double u3 = uniform(engine);
                    double u4 = uniform(engine);

                    double R3 = sqrt(-2 * log(u3));

                    double dx2 = R3 * cos(2 * M_PI * u3) * ouverture;
                    double dy2 = R3 * sin(2 * M_PI * u4) * ouverture;

                    Vector cameraPrime = camera + Vector(dx2, dy2, 0.0);
                    Vector uprime = camera + distanceMiseAuPoint * u - cameraPrime;
                    uprime.normalize();

                    rayonCamera = Rayon(cameraPrime, uprime);
                }

                lumiere += scene.getColor(rayonCamera, nRebonds);
            }
            lumiere = lumiere / numberRays;

            image[(i * W + j) * 3 + 0] = std::min(255., std::pow(lumiere[0], 0.45)); //  RED
            image[(i * W + j) * 3 + 1] = std::min(255., std::pow(lumiere[1], 0.45)); // GREEN
            image[(i * W + j) * 3 + 2] = std::min(255., std::pow(lumiere[2], 0.45)); // BLUE
        }
    }
    stbi_write_png("./images/damier.png", W, H, 3, &image[0], 0);
    return 0;
}
