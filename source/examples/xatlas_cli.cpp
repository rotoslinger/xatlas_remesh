/*
MIT License

Copyright (c) 2018-2020 Jonathan Young

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
*/
/*
xatlas_cli - Command Line Interface for xatlas UV Atlas Generator

This tool provides a comprehensive command line interface for the xatlas library
to generate UV atlases from 3D models with full control over all settings.

Input: an .obj model file.

Output:
	* an .obj model file with new UV coordinates
	* texture coordinates rasterized to images, colored by chart and by triangle
*/
#include <mutex>
#include <assert.h>
#include <stdarg.h>
#include <stdio.h>
#include <time.h>
#include <string>
#include <vector>
#include <map>

#include <stb_image_write.h>

#ifdef _MSC_VER
#pragma warning(push)
#pragma warning(disable : 4201)
#endif
#include <tiny_obj_loader.h>
#ifdef _MSC_VER
#pragma warning(pop)
#endif

#include <xatlas.h>

#ifdef _MSC_VER
#define FOPEN(_file, _filename, _mode) { if (fopen_s(&_file, _filename, _mode) != 0) _file = NULL; }
#define STRICMP _stricmp
#else
#define FOPEN(_file, _filename, _mode) _file = fopen(_filename, _mode)
#include <strings.h>
#define STRICMP strcasecmp
#endif

#define OBJ_TRIANGULATE 1 // Pass tinyobj::triangulation flag to tinyobjloader and treat all geometry as triangles.

static bool s_verbose = false;

// Command line options structure
struct CommandLineOptions {
    // Chart options
    float maxChartArea = 0.0f;
    float maxBoundaryLength = 0.0f;
    float normalDeviationWeight = 2.0f;
    float roundnessWeight = 0.01f;
    float straightnessWeight = 6.0f;
    float normalSeamWeight = 4.0f;
    float textureSeamWeight = 0.5f;
    float maxCost = 2.0f;
    uint32_t maxIterations = 1;
    bool useInputMeshUvs = false;
    bool fixWinding = false;
    
    // Pack options
    uint32_t maxChartSize = 0;
    uint32_t padding = 0;
    float texelsPerUnit = 0.0f;
    uint32_t resolution = 0;
    bool bilinear = true;
    bool blockAlign = false;
    bool bruteForce = false;
    bool createImage = false;
    bool rotateChartsToAxis = true;
    bool rotateCharts = true;
    
    // Output options
    std::string outputObj = "output.obj";
    std::string outputPrefix = "atlas";
    bool help = false;
};

class Stopwatch
{
public:
	Stopwatch() { reset(); }
	void reset() { m_start = clock(); }
	double elapsed() const { return (clock() - m_start) * 1000.0 / CLOCKS_PER_SEC; }
private:
	clock_t m_start;
};

static int Print(const char *format, ...)
{
	va_list arg;
	va_start(arg, format);
	printf("\r"); // Clear progress text.
	const int result = vprintf(format, arg);
	va_end(arg);
	return result;
}

// May be called from any thread.
static bool ProgressCallback(xatlas::ProgressCategory category, int progress, void *userData)
{
	// Don't interupt verbose printing.
	if (s_verbose)
		return true;
	Stopwatch *stopwatch = (Stopwatch *)userData;
	static std::mutex progressMutex;
	std::unique_lock<std::mutex> lock(progressMutex);
	if (progress == 0)
		stopwatch->reset();
	printf("\r   %s [", xatlas::StringForEnum(category));
	for (int i = 0; i < 10; i++)
		printf(progress / ((i + 1) * 10) ? "*" : " ");
	printf("] %d%%", progress);
	fflush(stdout);
	if (progress == 100)
		printf("\n      %.2f seconds (%g ms) elapsed\n", stopwatch->elapsed() / 1000.0, stopwatch->elapsed());
	return true;
}

static void RandomColor(uint8_t *color)
{
	for (int i = 0; i < 3; i++)
		color[i] = uint8_t((rand() % 255 + 192) * 0.5f);
}

static void SetPixel(uint8_t *dest, int destWidth, int x, int y, const uint8_t *color)
{
	uint8_t *pixel = &dest[x * 3 + y * (destWidth * 3)];
	pixel[0] = color[0];
	pixel[1] = color[1];
	pixel[2] = color[2];
}

// https://github.com/miloyip/line/blob/master/line_bresenham.c
// License: public domain.
static void RasterizeLine(uint8_t *dest, int destWidth, const int *p1, const int *p2, const uint8_t *color)
{
	const int dx = abs(p2[0] - p1[0]), sx = p1[0] < p2[0] ? 1 : -1;
	const int dy = abs(p2[1] - p1[1]), sy = p1[1] < p2[1] ? 1 : -1;
	int err = (dx > dy ? dx : -dy) / 2;
	int current[2];
	current[0] = p1[0];
	current[1] = p1[1];
	while (SetPixel(dest, destWidth, current[0], current[1], color), current[0] != p2[0] || current[1] != p2[1])
	{
		const int e2 = err;
		if (e2 > -dx) { err -= dy; current[0] += sx; }
		if (e2 < dy) { err += dx; current[1] += sy; }
	}
}

/*
https://github.com/ssloy/tinyrenderer/wiki/Lesson-2:-Triangle-rasterization-and-back-face-culling
Copyright Dmitry V. Sokolov

This software is provided 'as-is', without any express or implied warranty.
In no event will the authors be held liable for any damages arising from the use of this software.
Permission is granted to anyone to use this software for any purpose,
including commercial applications, and to alter it and redistribute it freely,
subject to the following restrictions:

1. The origin of this software must not be misrepresented; you must not claim that you wrote the original software. If you use this software in a product, an acknowledgment in the product documentation would be appreciated but is not required.
2. Altered source versions must be plainly marked as such, and must not be misrepresented as being the original software.
3. This notice may not be removed or altered from any source distribution.
*/
static void RasterizeTriangle(uint8_t *dest, int destWidth, const int *t0, const int *t1, const int *t2, const uint8_t *color)
{
	if (t0[1] > t1[1]) std::swap(t0, t1);
	if (t0[1] > t2[1]) std::swap(t0, t2);
	if (t1[1] > t2[1]) std::swap(t1, t2);
	int total_height = t2[1] - t0[1];
	for (int i = 0; i < total_height; i++) {
		bool second_half = i > t1[1] - t0[1] || t1[1] == t0[1];
		int segment_height = second_half ? t2[1] - t1[1] : t1[1] - t0[1];
		float alpha = (float)i / total_height;
		float beta = (float)(i - (second_half ? t1[1] - t0[1] : 0)) / segment_height;
		int A[2], B[2];
		for (int j = 0; j < 2; j++) {
			A[j] = int(t0[j] + (t2[j] - t0[j]) * alpha);
			B[j] = int(second_half ? t1[j] + (t2[j] - t1[j]) * beta : t0[j] + (t1[j] - t0[j]) * beta);
		}
		if (A[0] > B[0]) std::swap(A, B);
		for (int j = A[0]; j <= B[0]; j++)
			SetPixel(dest, destWidth, j, t0[1] + i, color);
	}
}



void PrintHelp(const char* programName)
{
    printf("xatlas_cli - Command Line Interface for xatlas UV Atlas Generator\n");
    printf("Usage: %s input_file.obj [options]\n\n", programName);
    printf("INPUT:\n");
    printf("  input_file.obj    Input 3D model file (.obj format)\n\n");
    
    printf("CHART OPTIONS:\n");
    printf("  --max-chart-area <float>        Max chart area (0 = no limit) [default: 0.0]\n");
    printf("  --max-boundary-length <float>   Max chart boundary length (0 = no limit) [default: 0.0]\n");
    printf("  --normal-deviation-weight <float> Weight for normal deviation [default: 2.0]\n");
    printf("  --roundness-weight <float>      Weight for chart roundness [default: 0.01]\n");
    printf("  --straightness-weight <float>   Weight for chart straightness [default: 6.0]\n");
    printf("  --normal-seam-weight <float>    Weight for normal seams [default: 4.0]\n");
    printf("  --texture-seam-weight <float>   Weight for texture seams [default: 0.5]\n");
    printf("  --max-cost <float>              Max cost for chart growth [default: 2.0]\n");
    printf("  --max-iterations <int>          Max chart growing iterations [default: 1]\n");
    printf("  --use-input-uvs                 Use input mesh UVs for charts\n");
    printf("  --fix-winding                   Fix texture coordinate winding\n\n");
    
    printf("PACK OPTIONS:\n");
    printf("  --max-chart-size <int>          Max chart size in pixels (0 = no limit) [default: 0]\n");
    printf("  --padding <int>                 Padding between charts in pixels [default: 0]\n");
    printf("  --texels-per-unit <float>       Texels per unit (0 = auto) [default: 0.0]\n");
    printf("  --resolution <int>              Atlas resolution (0 = auto) [default: 0]\n");
    printf("  --no-bilinear                   Disable bilinear filtering space\n");
    printf("  --block-align                   Align charts to 4x4 blocks\n");
    printf("  --brute-force                   Use brute force packing (slower but better)\n");
    printf("  --no-rotate-charts              Disable chart rotation\n");
    printf("  --no-rotate-to-axis             Disable rotation to convex hull axis\n\n");
    
    printf("OUTPUT OPTIONS:\n");
    printf("  --output-obj <filename>         Output OBJ filename [default: output.obj]\n");
    printf("  --output-prefix <prefix>        Output image prefix [default: atlas]\n");
    printf("  --verbose                       Enable verbose output\n");
    printf("  --help, -h                      Show this help message\n\n");
    
    printf("EXAMPLES:\n");
    printf("  %s model.obj\n", programName);
    printf("  %s model.obj --padding 4 --resolution 2048 --verbose\n", programName);
    printf("  %s model.obj --max-chart-area 100 --max-cost 1.5 --texels-per-unit 64\n", programName);
    printf("  %s model.obj --output-obj my_output.obj --output-prefix my_atlas\n", programName);
}

bool ParseCommandLine(int argc, char* argv[], CommandLineOptions& options, std::string& inputFile)
{
    if (argc < 2) {
        PrintHelp(argv[0]);
        return false;
    }
    
    // Check for help flag first
    for (int i = 1; i < argc; i++) {
        std::string arg = argv[i];
        if (arg == "--help" || arg == "-h") {
            options.help = true;
            return true;
        }
    }
    
    inputFile = argv[1];
    
    for (int i = 2; i < argc; i++) {
        std::string arg = argv[i];
        
        if (arg == "--help" || arg == "-h") {
            options.help = true;
            return true;
        }
        else if (arg == "--verbose") {
            s_verbose = true;
        }
        else if (arg == "--max-chart-area" && i + 1 < argc) {
            options.maxChartArea = std::stof(argv[++i]);
        }
        else if (arg == "--max-boundary-length" && i + 1 < argc) {
            options.maxBoundaryLength = std::stof(argv[++i]);
        }
        else if (arg == "--normal-deviation-weight" && i + 1 < argc) {
            options.normalDeviationWeight = std::stof(argv[++i]);
        }
        else if (arg == "--roundness-weight" && i + 1 < argc) {
            options.roundnessWeight = std::stof(argv[++i]);
        }
        else if (arg == "--straightness-weight" && i + 1 < argc) {
            options.straightnessWeight = std::stof(argv[++i]);
        }
        else if (arg == "--normal-seam-weight" && i + 1 < argc) {
            options.normalSeamWeight = std::stof(argv[++i]);
        }
        else if (arg == "--texture-seam-weight" && i + 1 < argc) {
            options.textureSeamWeight = std::stof(argv[++i]);
        }
        else if (arg == "--max-cost" && i + 1 < argc) {
            options.maxCost = std::stof(argv[++i]);
        }
        else if (arg == "--max-iterations" && i + 1 < argc) {
            options.maxIterations = std::stoul(argv[++i]);
        }
        else if (arg == "--use-input-uvs") {
            options.useInputMeshUvs = true;
        }
        else if (arg == "--fix-winding") {
            options.fixWinding = true;
        }
        else if (arg == "--max-chart-size" && i + 1 < argc) {
            options.maxChartSize = std::stoul(argv[++i]);
        }
        else if (arg == "--padding" && i + 1 < argc) {
            options.padding = std::stoul(argv[++i]);
        }
        else if (arg == "--texels-per-unit" && i + 1 < argc) {
            options.texelsPerUnit = std::stof(argv[++i]);
        }
        else if (arg == "--resolution" && i + 1 < argc) {
            options.resolution = std::stoul(argv[++i]);
        }
        else if (arg == "--no-bilinear") {
            options.bilinear = false;
        }
        else if (arg == "--block-align") {
            options.blockAlign = true;
        }
        else if (arg == "--brute-force") {
            options.bruteForce = true;
        }
        else if (arg == "--no-rotate-charts") {
            options.rotateCharts = false;
        }
        else if (arg == "--no-rotate-to-axis") {
            options.rotateChartsToAxis = false;
        }
        else if (arg == "--output-obj" && i + 1 < argc) {
            options.outputObj = argv[++i];
        }
        else if (arg == "--output-prefix" && i + 1 < argc) {
            options.outputPrefix = argv[++i];
        }
        else {
            printf("Error: Unknown option '%s'\n", arg.c_str());
            printf("Use --help for usage information\n");
            return false;
        }
    }
    
    return true;
}

int main(int argc, char *argv[])
{
    CommandLineOptions options;
    std::string inputFile;
    
    if (!ParseCommandLine(argc, argv, options, inputFile)) {
        return EXIT_FAILURE;
    }
    
    if (options.help) {
        PrintHelp(argv[0]);
        return EXIT_SUCCESS;
    }
    
    // Load object file.
    printf("Loading '%s'...\n", inputFile.c_str());
    std::vector<tinyobj::shape_t> shapes;
    std::vector<tinyobj::material_t> materials;
    std::string err;
    unsigned int flags = 0;
#if OBJ_TRIANGULATE
    flags = tinyobj::triangulation;
#endif
    if (!tinyobj::LoadObj(shapes, materials, err, inputFile.c_str(), NULL, flags)) {
        printf("Error: %s\n", err.c_str());
        return EXIT_FAILURE;
    }
    if (shapes.size() == 0) {
        printf("Error: no shapes in obj file\n");
        return EXIT_FAILURE;
    }
    printf("   %d shapes\n", (int)shapes.size());
    
    // Create empty atlas.
    xatlas::SetPrint(Print, s_verbose);
    xatlas::Atlas *atlas = xatlas::Create();
    
    // Set progress callback.
    Stopwatch globalStopwatch, stopwatch;
    xatlas::SetProgressCallback(atlas, ProgressCallback, &stopwatch);
    
    // Add meshes to atlas.
    uint32_t totalVertices = 0, totalFaces = 0;
    for (int i = 0; i < (int)shapes.size(); i++) {
        const tinyobj::mesh_t &objMesh = shapes[i].mesh;
        xatlas::MeshDecl meshDecl;
        meshDecl.vertexCount = (uint32_t)objMesh.positions.size() / 3;
        meshDecl.vertexPositionData = objMesh.positions.data();
        meshDecl.vertexPositionStride = sizeof(float) * 3;
        if (!objMesh.normals.empty()) {
            meshDecl.vertexNormalData = objMesh.normals.data();
            meshDecl.vertexNormalStride = sizeof(float) * 3;
        }
        if (!objMesh.texcoords.empty()) {
            meshDecl.vertexUvData = objMesh.texcoords.data();
            meshDecl.vertexUvStride = sizeof(float) * 2;
        }
        meshDecl.indexCount = (uint32_t)objMesh.indices.size();
        meshDecl.indexData = objMesh.indices.data();
        meshDecl.indexFormat = xatlas::IndexFormat::UInt32;

        xatlas::AddMeshError error = xatlas::AddMesh(atlas, meshDecl, (uint32_t)shapes.size());
        if (error != xatlas::AddMeshError::Success) {
            xatlas::Destroy(atlas);
            printf("\rError adding mesh %d '%s': %s\n", i, shapes[i].name.c_str(), xatlas::StringForEnum(error));
            return EXIT_FAILURE;
        }
        totalVertices += meshDecl.vertexCount;
        if (meshDecl.faceCount > 0)
            totalFaces += meshDecl.faceCount;
        else
            totalFaces += meshDecl.indexCount / 3; // Assume triangles if MeshDecl::faceCount not specified.
    }
    xatlas::AddMeshJoin(atlas); // Not necessary. Only called here so geometry totals are printed after the AddMesh progress indicator.
    printf("   %u total vertices\n", totalVertices);
    printf("   %u total faces\n", totalFaces);
    
    // Configure chart options
    xatlas::ChartOptions chartOptions;
    chartOptions.maxChartArea = options.maxChartArea;
    chartOptions.maxBoundaryLength = options.maxBoundaryLength;
    chartOptions.normalDeviationWeight = options.normalDeviationWeight;
    chartOptions.roundnessWeight = options.roundnessWeight;
    chartOptions.straightnessWeight = options.straightnessWeight;
    chartOptions.normalSeamWeight = options.normalSeamWeight;
    chartOptions.textureSeamWeight = options.textureSeamWeight;
    chartOptions.maxCost = options.maxCost;
    chartOptions.maxIterations = options.maxIterations;
    chartOptions.useInputMeshUvs = options.useInputMeshUvs;
    chartOptions.fixWinding = options.fixWinding;
    
    // Configure pack options
    xatlas::PackOptions packOptions;
    packOptions.maxChartSize = options.maxChartSize;
    packOptions.padding = options.padding;
    packOptions.texelsPerUnit = options.texelsPerUnit;
    packOptions.resolution = options.resolution;
    packOptions.bilinear = options.bilinear;
    packOptions.blockAlign = options.blockAlign;
    packOptions.bruteForce = options.bruteForce;
    packOptions.createImage = options.createImage;
    packOptions.rotateChartsToAxis = options.rotateChartsToAxis;
    packOptions.rotateCharts = options.rotateCharts;
    
    // Generate atlas with custom options
    printf("Generating atlas with custom options...\n");
    if (s_verbose) {
        printf("Chart options: maxArea=%.1f, maxBoundary=%.1f, normalDev=%.1f, roundness=%.3f, straightness=%.1f, normalSeam=%.1f, textureSeam=%.1f, maxCost=%.1f, iterations=%u\n",
               chartOptions.maxChartArea, chartOptions.maxBoundaryLength, chartOptions.normalDeviationWeight,
               chartOptions.roundnessWeight, chartOptions.straightnessWeight, chartOptions.normalSeamWeight,
               chartOptions.textureSeamWeight, chartOptions.maxCost, chartOptions.maxIterations);
        printf("Pack options: maxSize=%u, padding=%u, texelsPerUnit=%.1f, resolution=%u, bilinear=%s, blockAlign=%s, bruteForce=%s, rotateCharts=%s, rotateToAxis=%s\n",
               packOptions.maxChartSize, packOptions.padding, packOptions.texelsPerUnit, packOptions.resolution,
               packOptions.bilinear ? "true" : "false", packOptions.blockAlign ? "true" : "false",
               packOptions.bruteForce ? "true" : "false", packOptions.rotateCharts ? "true" : "false",
               packOptions.rotateChartsToAxis ? "true" : "false");
    }
    
    xatlas::Generate(atlas, chartOptions, packOptions);
    printf("   %d charts\n", atlas->chartCount);
    printf("   %d atlases\n", atlas->atlasCount);
    for (uint32_t i = 0; i < atlas->atlasCount; i++)
        printf("      %d: %0.2f%% utilization\n", i, atlas->utilization[i] * 100.0f);
    printf("   %ux%u resolution\n", atlas->width, atlas->height);
    totalVertices = 0;
    for (uint32_t i = 0; i < atlas->meshCount; i++) {
        const xatlas::Mesh &mesh = atlas->meshes[i];
        totalVertices += mesh.vertexCount;
        // Input and output index counts always match.
        assert(mesh.indexCount == (uint32_t)shapes[i].mesh.indices.size());
    }
    printf("   %u total vertices\n", totalVertices);
    printf("%.2f seconds (%g ms) elapsed total\n", globalStopwatch.elapsed() / 1000.0, globalStopwatch.elapsed());
    
    // Write meshes.
    printf("Writing '%s'...\n", options.outputObj.c_str());
    FILE *file;
    FOPEN(file, options.outputObj.c_str(), "w");
    if (file) {
        uint32_t firstVertex = 0;
        for (uint32_t i = 0; i < atlas->meshCount; i++) {
            const xatlas::Mesh &mesh = atlas->meshes[i];
            for (uint32_t v = 0; v < mesh.vertexCount; v++) {
                const xatlas::Vertex &vertex = mesh.vertexArray[v];
                const float *pos = &shapes[i].mesh.positions[vertex.xref * 3];
                fprintf(file, "v %g %g %g\n", pos[0], pos[1], pos[2]);
                if (!shapes[i].mesh.normals.empty()) {
                    const float *normal = &shapes[i].mesh.normals[vertex.xref * 3];
                    fprintf(file, "vn %g %g %g\n", normal[0], normal[1], normal[2]);
                }
                fprintf(file, "vt %g %g\n", vertex.uv[0] / atlas->width, vertex.uv[1] / atlas->height);
            }
            fprintf(file, "o %s\n", shapes[i].name.c_str());
            fprintf(file, "s off\n");
            for (uint32_t f = 0; f < mesh.indexCount; f += 3) {
                fprintf(file, "f ");
                for (uint32_t j = 0; j < 3; j++) {
                    const uint32_t index = firstVertex + mesh.indexArray[f + j] + 1; // 1-indexed
                    fprintf(file, "%d/%d/%d%c", index, index, index, j == 2 ? '\n' : ' ');
                }
            }
            firstVertex += mesh.vertexCount;
        }
        fclose(file);
    }
    
    if (atlas->width > 0 && atlas->height > 0) {
        printf("Rasterizing result...\n");
        // Dump images.
        std::vector<uint8_t> outputTrisImage, outputChartsImage;
        const uint32_t imageDataSize = atlas->width * atlas->height * 3;
        outputTrisImage.resize(atlas->atlasCount * imageDataSize);
        outputChartsImage.resize(atlas->atlasCount * imageDataSize);
        for (uint32_t i = 0; i < atlas->meshCount; i++) {
            const xatlas::Mesh &mesh = atlas->meshes[i];
            // Rasterize mesh triangles.
            const uint8_t white[] = { 255, 255, 255 };
            const uint32_t faceCount = mesh.indexCount / 3;
            uint32_t faceFirstIndex = 0;
            for (uint32_t f = 0; f < faceCount; f++) {
                int32_t atlasIndex = -1;
                int verts[255][2];
                const uint32_t faceVertexCount = 3;
                for (uint32_t j = 0; j < faceVertexCount; j++) {
                    const xatlas::Vertex &v = mesh.vertexArray[mesh.indexArray[faceFirstIndex + j]];
                    atlasIndex = v.atlasIndex; // The same for every vertex in the face.
                    verts[j][0] = int(v.uv[0]);
                    verts[j][1] = int(v.uv[1]);
                }
                if (atlasIndex < 0)
                    continue; // Skip faces that weren't atlased.
                uint8_t color[3];
                RandomColor(color);
                uint8_t *imageData = &outputTrisImage[atlasIndex * imageDataSize];
                RasterizeTriangle(imageData, atlas->width, verts[0], verts[1], verts[2], color);
                for (uint32_t j = 0; j < faceVertexCount; j++)
                    RasterizeLine(imageData, atlas->width, verts[j], verts[(j + 1) % faceVertexCount], white);
                faceFirstIndex += faceVertexCount;
            }
            // Rasterize mesh charts.
            for (uint32_t j = 0; j < mesh.chartCount; j++) {
                const xatlas::Chart *chart = &mesh.chartArray[j];
                uint8_t color[3];
                RandomColor(color);
                for (uint32_t k = 0; k < chart->faceCount; k++) {
                    const uint32_t face = chart->faceArray[k];
                    const uint32_t faceVertexCount = 3;
                    faceFirstIndex = face * 3;
                    int verts[255][2];
                    for (uint32_t l = 0; l < faceVertexCount; l++) {
                        const xatlas::Vertex &v = mesh.vertexArray[mesh.indexArray[faceFirstIndex + l]];
                        verts[l][0] = int(v.uv[0]);
                        verts[l][1] = int(v.uv[1]);
                    }
                    uint8_t *imageData = &outputChartsImage[chart->atlasIndex * imageDataSize];
                    RasterizeTriangle(imageData, atlas->width, verts[0], verts[1], verts[2], color);
                    for (uint32_t l = 0; l < faceVertexCount; l++)
                        RasterizeLine(imageData, atlas->width, verts[l], verts[(l + 1) % faceVertexCount], white);
                }
            }
        }
        for (uint32_t i = 0; i < atlas->atlasCount; i++) {
            char filename[256];
            snprintf(filename, sizeof(filename), "%s_tris%02u.tga", options.outputPrefix.c_str(), i);
            printf("Writing '%s'...\n", filename);
            stbi_write_tga(filename, atlas->width, atlas->height, 3, &outputTrisImage[i * imageDataSize]);
            snprintf(filename, sizeof(filename), "%s_charts%02u.tga", options.outputPrefix.c_str(), i);
            printf("Writing '%s'...\n", filename);
            stbi_write_tga(filename, atlas->width,atlas->height, 3, &outputChartsImage[i * imageDataSize]);
        }
    }
    // Cleanup.
    xatlas::Destroy(atlas);
    printf("Done\n");
    return EXIT_SUCCESS;
}
