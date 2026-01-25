/* tslint:disable */
/* eslint-disable */

export class BamFile {
    free(): void;
    [Symbol.dispose](): void;
    calculate_coverage(ref_name: string, start: number, end: number): any;
    generate_pileup(ref_name: string, start: number, end: number, min_mapq: number, min_baseq: number): any;
    get_read_count(): bigint;
    get_reference_names(): string[];
    get_region_stats(ref_name: string, start: number, end: number): any;
    constructor(data: Uint8Array);
}

export class ChromatiqEngine {
    free(): void;
    [Symbol.dispose](): void;
    calculate_combined_coverage(ref_name: string, start: number, end: number): any;
    get_all_variants(): any;
    get_bam_file_count(): number;
    get_vcf_file_count(): number;
    load_bam(data: Uint8Array): number;
    load_vcf(data: Uint8Array): number;
    constructor();
    unload_bam(index: number): void;
    unload_vcf(index: number): void;
}

export class VcfFile {
    free(): void;
    [Symbol.dispose](): void;
    get_high_quality_variants(min_qual: number): any;
    get_passed_variants(): any;
    get_stats(): any;
    get_variants(): any;
    get_variants_in_region(chrom: string, start: number, end: number): any;
    constructor(data: Uint8Array);
}

export function get_capabilities(): string;

export function get_info(): string;

export function get_version(): string;

export function init(): void;

export type InitInput = RequestInfo | URL | Response | BufferSource | WebAssembly.Module;

export interface InitOutput {
    readonly memory: WebAssembly.Memory;
    readonly __wbg_bamfile_free: (a: number, b: number) => void;
    readonly __wbg_chromatiqengine_free: (a: number, b: number) => void;
    readonly __wbg_vcffile_free: (a: number, b: number) => void;
    readonly bamfile_calculate_coverage: (a: number, b: number, c: number, d: number, e: number) => [number, number, number];
    readonly bamfile_generate_pileup: (a: number, b: number, c: number, d: number, e: number, f: number, g: number) => [number, number, number];
    readonly bamfile_get_read_count: (a: number) => bigint;
    readonly bamfile_get_reference_names: (a: number) => [number, number];
    readonly bamfile_new: (a: number, b: number) => [number, number, number];
    readonly chromatiqengine_calculate_combined_coverage: (a: number, b: number, c: number, d: number, e: number) => [number, number, number];
    readonly chromatiqengine_get_all_variants: (a: number) => [number, number, number];
    readonly chromatiqengine_get_bam_file_count: (a: number) => number;
    readonly chromatiqengine_get_vcf_file_count: (a: number) => number;
    readonly chromatiqengine_load_bam: (a: number, b: number, c: number) => [number, number, number];
    readonly chromatiqengine_load_vcf: (a: number, b: number, c: number) => [number, number, number];
    readonly chromatiqengine_new: () => number;
    readonly chromatiqengine_unload_bam: (a: number, b: number) => [number, number];
    readonly chromatiqengine_unload_vcf: (a: number, b: number) => [number, number];
    readonly get_capabilities: () => [number, number];
    readonly get_info: () => [number, number];
    readonly get_version: () => [number, number];
    readonly vcffile_get_high_quality_variants: (a: number, b: number) => [number, number, number];
    readonly vcffile_get_passed_variants: (a: number) => [number, number, number];
    readonly vcffile_get_stats: (a: number) => [number, number, number];
    readonly vcffile_get_variants: (a: number) => [number, number, number];
    readonly vcffile_get_variants_in_region: (a: number, b: number, c: number, d: number, e: number) => [number, number, number];
    readonly vcffile_new: (a: number, b: number) => [number, number, number];
    readonly init: () => void;
    readonly bamfile_get_region_stats: (a: number, b: number, c: number, d: number, e: number) => [number, number, number];
    readonly __wbindgen_malloc: (a: number, b: number) => number;
    readonly __wbindgen_realloc: (a: number, b: number, c: number, d: number) => number;
    readonly __wbindgen_exn_store: (a: number) => void;
    readonly __externref_table_alloc: () => number;
    readonly __wbindgen_externrefs: WebAssembly.Table;
    readonly __wbindgen_free: (a: number, b: number, c: number) => void;
    readonly __externref_table_dealloc: (a: number) => void;
    readonly __externref_drop_slice: (a: number, b: number) => void;
    readonly __wbindgen_start: () => void;
}

export type SyncInitInput = BufferSource | WebAssembly.Module;

/**
 * Instantiates the given `module`, which can either be bytes or
 * a precompiled `WebAssembly.Module`.
 *
 * @param {{ module: SyncInitInput }} module - Passing `SyncInitInput` directly is deprecated.
 *
 * @returns {InitOutput}
 */
export function initSync(module: { module: SyncInitInput } | SyncInitInput): InitOutput;

/**
 * If `module_or_path` is {RequestInfo} or {URL}, makes a request and
 * for everything else, calls `WebAssembly.instantiate` directly.
 *
 * @param {{ module_or_path: InitInput | Promise<InitInput> }} module_or_path - Passing `InitInput` directly is deprecated.
 *
 * @returns {Promise<InitOutput>}
 */
export default function __wbg_init (module_or_path?: { module_or_path: InitInput | Promise<InitInput> } | InitInput | Promise<InitInput>): Promise<InitOutput>;
