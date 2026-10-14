/******************************************************************************************************
// mnovaBatchExportJson.qs
// 用途: 无界面模式把命令行传入的 .mnova 文档中每条 NMR 谱导出为 JSON, 供 NMR2Structure 主流程桥接读取.
// 调用: MestReNova.exe --nogui -w <输入.mnova...> --sf mnovaBatchExportJson
// 输出: <输入.mnova 所在目录>/mnova_export/<文档名>/<文档名>__page<NN>__item<K>__<核素[_2D]>.json
//       (输出根从文档自身路径动态推导, 不在脚本中写死任何路径)
//       JSON 字段: title, dim, nucleus(/nucleus2), frequency(/2)_MHz, solvent, firstXppm/lastXppm,
//                  pointCount, peaks[{ppm1,(ppm2),intensity,area,annotation}], data(仅1D, real 数组)
/*****************************************************************************************************/
/*globals mainWindow, serialization, NMRSpectrum, File, TextStream, Dir, Document, Application, print, JSON*/

function mnovaBatchExportJson() {
    "use strict";

    var docs, di, doc, pageCount, pageIndex, itemIndex, page, itemCount, item, spec,
        outRoot, outDir, docName, docPath, exported = 0, summary = [];

    function basename(path) {
        return path.substring(Math.max(path.lastIndexOf("/"), path.lastIndexOf("\\")) + 1)
            .replace(/\.[^.]+$/, "");
    }

    // 文档完整路径去掉文件名 → 所在目录 (doc.name 返回 "E:/.../clxh98.mnova" 全路径)
    function dirname(path) {
        return path.substring(0, Math.max(path.lastIndexOf("/"), path.lastIndexOf("\\")));
    }

    function nucleusLabel(spec3, dim) {
        var nuc = spec3.nucleus(dim, false);
        return String(nuc).replace(/[^0-9A-Za-z]/g, "");
    }

    function expName(spec3) {
        var n1, d = spec3.dimCount;
        n1 = nucleusLabel(spec3, 1);
        if (d === 1) { return n1; }
        return n1 + "_" + nucleusLabel(spec3, 2) + "_" + d + "D";
    }

    function exportSpectrum(spec3, outDir3, base3, pageNo, itemNo) {
        var file3, strm3, obj3, peaks3, peak3, j, realArr3, p, objKey,
            nmrProc3;

        if (!spec3.isValid()) { return false; }

        // 峰表为空时自动拾峰: 1D 用 GSD, 2D 用 Standard (参照官方 DBPlugin.qs getNMRPeaks 做法)
        peaks3 = spec3.peaks();
        if (!peaks3 || peaks3.count === 0) {
            try {
                nmrProc3 = new NMRProcessing(spec3.proc);
                nmrProc3.setParameter("PP.Method", spec3.dimCount === 1 ? "GSD" : "Standard");
                nmrProc3.setParameter("PP.AutoClassify", true);
                nmrProc3.setParameter("PP.Apply", true);
                spec3.processCurSpec(nmrProc3);
                spec3.update();
                peaks3 = spec3.peaks();
            } catch (ePeak) {
                print("[WARN] autopick failed page" + pageNo + ": " + ePeak);
                peaks3 = spec3.peaks();
            }
        }

        obj3 = {};
        obj3.title = spec3.title;
        obj3.experimentType = spec3.experimentType;  // MNova 内部实验类型 (HSQC/HMBC/COSY/ROESY/DEPT-135/1H/13C...)
        obj3.originalFormat = spec3.originalFormat;
        obj3.dim = spec3.dimCount;
        obj3.solvent = spec3.solvent;
        obj3.nucleus = nucleusLabel(spec3, 1);
        obj3.frequency_MHz = spec3.frequency();
        if (spec3.dimCount >= 2) {
            obj3.nucleus2 = nucleusLabel(spec3, 2);
            obj3.frequency2_MHz = spec3.frequency(2);
        }
        try {
            obj3.firstXppm = spec3.unitsToUnits("pt", "ppm", 0, 1);
            obj3.lastXppm = spec3.unitsToUnits("pt", "ppm", spec3.count(1) - 1, 1);
        } catch (e1) {
            obj3.firstXppm = null;
            obj3.lastXppm = null;
        }
        obj3.pointCount = spec3.count(1);

        peaks3 = spec3.peaks();
        obj3.peaks = [];
        if (peaks3 && peaks3.count) {
            for (j = 0; j < peaks3.count; j++) {
                peak3 = peaks3.at(j);
                obj3.peaks.push({
                    ppm1: peak3.delta(1),
                    intensity: peak3.intensity,
                    area: peak3.integralInPts(spec3),
                    annotation: peak3.annotation
                });
                if (spec3.dimCount === 2) {
                    obj3.peaks[obj3.peaks.length - 1].ppm2 = peak3.delta(2);
                }
            }
        }

        // 1D 导出全实数矩阵; 2D 导出完整实数矩阵 (复核模块需在矩阵上找局部极大值, 如融合椭圆拆分)
        if (spec3.dimCount === 1) {
            realArr3 = [];
            for (p = 0; p < spec3.count(1); p++) {
                realArr3.push(spec3.real(p));
            }
            obj3.data = realArr3;
        } else if (spec3.dimCount === 2) {
            // 行=F1 (碳维, count(1)), 列=F2 (氢维, count(2)); real(i,j) 第 i 行第 j 列
            var nF13, nF23, i3, rowArr3;
            nF13 = spec3.count(1);
            nF23 = spec3.count(2);
            try {
                obj3.firstYppm = spec3.unitsToUnits("pt", "ppm", 0, 2);
                obj3.lastYppm = spec3.unitsToUnits("pt", "ppm", nF13 - 1, 2);
            } catch (e2d) {
                obj3.firstYppm = null;
                obj3.lastYppm = null;
            }
            obj3.pointCount2 = nF23;
            var mat3 = [];
            for (i3 = 0; i3 < nF13; i3++) {
                rowArr3 = [];
                for (p = 0; p < nF23; p++) {
                    rowArr3.push(spec3.real(i3, p));
                }
                mat3.push(rowArr3);
            }
            obj3.data2d = mat3;
        }

        objKey = base3 + "__page" + pageNo + "__item" + itemNo + "__" + expName(spec3);
        file3 = new File(outDir3 + "/" + objKey + ".json");
        if (!file3.open(File.WriteOnly)) { return false; }
        strm3 = new TextStream(file3);
        strm3.write(JSON.stringify(obj3));
        strm3.flush();
        file3.close();
        print("[OK] " + objKey + " peaks=" + obj3.peaks.length +
              (obj3.data ? " points=" + obj3.data.length : ""));
        return true;
    }

    docs = mainWindow.documents();
    for (di = 0; di < docs.length; di++) {
        doc = docs[di];
        // 跳过空白辅助文档
        if (!doc || doc.pageCount() === 0) { continue; }

        docName = basename(String(doc.name));
        docPath = String(doc.name);
        // 输出根: <文档所在目录>/mnova_export (动态推导自文档路径, 不写死)
        outRoot = dirname(docPath) + "/mnova_export";
        outDir = outRoot + "/" + docName;
        new Dir(outDir).mkpath(outDir);

        pageCount = doc.pageCount();
        print("[INFO] doc=" + docName + " pages=" + pageCount);
        for (pageIndex = 0; pageIndex < pageCount; pageIndex++) {
            page = doc.page(pageIndex);
            itemCount = page.itemCount();
            for (itemIndex = 0; itemIndex < itemCount; itemIndex++) {
                item = page.item(itemIndex);
                if (item.name === "NMR Spectrum") {
                    spec = new NMRSpectrum(item);
                    if (exportSpectrum(spec, outDir, docName, pageIndex + 1, itemIndex + 1)) {
                        exported++;
                    }
                }
            }
        }
        summary.push(docName + "=" + outDir);
    }
    print("[DONE] exported=" + exported);
    // 汇总落盘到输出根目录 (最后一个文档所在目录/mnova_export), 供外部脚本读取导出位置
    var sumFile = new File(outRoot + "/_summary.txt");
    if (sumFile.open(File.WriteOnly)) {
        var sumStrm = new TextStream(sumFile);
        for (var si = 0; si < summary.length; si++) {
            sumStrm.write(summary[si] + "\n");
        }
        sumStrm.flush();
        sumFile.close();
    }
    return "OK exported=" + exported;
}
