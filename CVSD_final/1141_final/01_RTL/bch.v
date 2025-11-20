module bch(
	input clk,
	input rstn,
	input mode,
	input [1:0] code,
	input set,
	input [63:0] idata,
	output ready,
	output finish,
	output [9:0] odata
);
// -------------------------------------------------------------------------------------------------------------------------- //
// Declaration
// -------------------------------------------------------------------------------------------------------------------------- //
	reg       ready_reg;
	assign    ready = ready_reg;

	reg [1:0] code_reg; // 1: (63,51); 2: (255,239); 3: (1023,983)
	reg       mode_reg; // 0: hard-decision decoding; 1: soft-decision decoding

// counter
	reg [2:0] cnt, cnt_n;

// LLR
	reg [1:1023] llr, llr_n; // 1: MSB ;  1023: LSB
	reg [7:0]    idata_in;

// CRC
	reg [0:4] w_CRC, r_CRC; // 0: MSB ; 4: LSB
	// reg [0:4] 
	reg [2:0] w_CRC_remainder_0, w_CRC_remainder_1, w_CRC_remainder_2, w_CRC_remainder_3, w_CRC_remainder_4, w_CRC_remainder_5, w_CRC_remainder_6, w_CRC_remainder_7, w_CRC_remainder_8, w_CRC_remainder_9, w_CRC_remainder_10;
    

// assign finish = (cnt == 7 || cnt == 8)? 1:0;

	reg [2:0] c_state, n_state;
	localparam S_IDLE  = 3'b000,
			   S_FETCH = 3'b001,
			   S_CRC   = 3'b011,
			   S_END   = 3'b111;


// -------------------------------------------------------------------------------------------------------------------------- //
// Logic
// -------------------------------------------------------------------------------------------------------------------------- //
// FSM
	always @(posedge clk) begin
		if(~rstn) begin
			c_state <= S_IDLE;
		end
		else begin
			c_state <= n_state;
		end
	end
	always @* begin
		n_state = c_state;
		case (c_state)
			S_IDLE: n_state = S_FETCH;
			S_FETCH: n_state = (cnt == 7)? S_CRC:S_FETCH;
			S_CRC  : n_state = (cnt == 7)? S_END:S_CRC;
			
			//default: 
		endcase
	end

// code info.
	always @(posedge clk) begin
		if(~rstn) begin
			code_reg <= 0;
			mode_reg <= 0;
		end
		else if(set) begin
			code_reg <= code;
			mode_reg <= mode;
		end
	end

// counter
	always @(posedge clk) begin
		if(~rstn) begin
			cnt <= 0;
		end
		else begin
			cnt <= cnt_n;
		end
	end

	always @* begin
		cnt_n = cnt;
		case(c_state)
			S_FETCH: cnt_n = (cnt == 7)? 0:cnt + 1;
		endcase
	end

// llr
	always @(posedge clk) begin
		if(~rstn) begin
			llr <= 0;
		end
		else begin
			llr <= llr_n;
		end
	end
	always @* begin
		llr_n = llr;
		case(c_state)
			S_FETCH: begin
				case(cnt)
					0: llr_n[1:7] = {idata[55], idata[47], idata[39], idata[31], idata[23], idata[15], idata[7]};
					1: llr_n[8:15] = {idata[63], idata[55], idata[47], idata[39], idata[31], idata[23], idata[15], idata[7]};
					2: llr_n[16:23] = {idata[63], idata[55], idata[47], idata[39], idata[31], idata[23], idata[15], idata[7]};
					3: llr_n[24:31] = {idata[63], idata[55], idata[47], idata[39], idata[31], idata[23], idata[15], idata[7]};
					4: llr_n[32:39] = {idata[63], idata[55], idata[47], idata[39], idata[31], idata[23], idata[15], idata[7]};
					5: llr_n[40:47] = {idata[63], idata[55], idata[47], idata[39], idata[31], idata[23], idata[15], idata[7]};
					6: llr_n[48:55] = {idata[63], idata[55], idata[47], idata[39], idata[31], idata[23], idata[15], idata[7]};
					7: llr_n[56:63] = {idata[63], idata[55], idata[47], idata[39], idata[31], idata[23], idata[15], idata[7]};
				endcase
			end
		endcase
	end









// output 
	always @* begin
		if(c_state == S_FETCH) ready_reg = 1;
		else ready_reg = 0;
	end

endmodule

